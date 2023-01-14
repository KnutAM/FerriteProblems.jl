using Ferrite, FerriteMeshParser, Tensors
using FerriteAssembly, FerriteProblems, FerriteNeumann, FESolvers
using MaterialModelsBase
import FerriteProblems as FP
import MaterialModelsBase as MMB

struct Elastic{T} <: AbstractMaterial
    G::T
    K::T
    E4::SymmetricTensor{4,2,T,9}
end
function Elastic(;E=2.e3, ν=0.3)
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    I2 = one(SymmetricTensor{2,2})
    I4vol = I2⊗I2
    I4dev = minorsymmetric(otimesu(I2,I2)) - I4vol / 3
    E4 = 2G*I4dev + K*I4vol
    return Elastic(G, K, E4)
end;

function MMB.material_response(m::Elastic, ϵ, args...; kwargs...)
    σ = m.E4 ⊡ ϵ
    return σ, m.E4, NoMaterialState()
end;

struct PoroElastic{E<:Elastic,T}
    elastic::E
    k::T    # [mm^4/Ns] Permeability
    α::T    # [-] Biot's coefficient
    β::T    # [1/MPa] Liquid bulk modulus
end
function PoroElastic(;elastic=Elastic(), k=0.05, α=1.0, β=1/2e3)
    return PoroElastic(elastic, k, α, β)
end;

function FerriteAssembly.element_residual!(re, state, ae, material::PoroElastic, cv::NamedTuple, dh_fh, Δt, buffer)
    # Setup cellvalues and give easier names
    cv_u = cv[:u]
    cv_p = cv[:p]
    num_u = getnbasefunctions(cv_u)
    num_p = getnbasefunctions(cv_p)

    # Assign views to the matrix and vector parts
    ae_old = FerriteAssembly.get_aeold(buffer)
    udofs = dof_range(dh_fh, :u)
    pdofs = dof_range(dh_fh, :p)
    ru = @view re[udofs]
    rp = @view re[pdofs]
    ue = @view ae[udofs]
    pe = @view ae[pdofs]
    ue_old = @view ae_old[udofs]
    pe_old = @view ae_old[pdofs]

    # Assemble stiffness and force vectors
    for q_point in 1:getnquadpoints(cv_u)
        # Calculate variables in the current quadrature point
        dΩ = getdetJdV(cv_u, q_point)
        ϵ = function_symmetric_gradient(cv_u, q_point, ue)
        ϵ_old = function_symmetric_gradient(cv_u, q_point, ue_old)
        p = function_value(cv_p, q_point, pe)
        p_old = function_value(cv_p, q_point, pe_old)
        ∇p = function_gradient(cv_p, q_point, pe)
        pdot = (p-p_old)/Δt
        div_udot = (tr(ϵ)-tr(ϵ_old))/Δt
        σeff = material.elastic.E4 ⊡ ϵ

        # Assemble residual contributions
        for iᵤ in 1:num_u
            ∇δNu = shape_symmetric_gradient(cv_u, q_point, iᵤ)
            div_δNu = shape_divergence(cv_u, q_point, iᵤ)
            ru[iᵤ] += (∇δNu ⊡ σeff - div_δNu*material.α*p)*dΩ
        end
        for iₚ in 1:num_p
            δNp = shape_value(cv_p, q_point, iₚ)
            ∇δNp = shape_gradient(cv_p, q_point, iₚ)
            rp[iₚ] += (δNp*(material.α*div_udot + material.β*pdot) + (∇δNp ⋅ ∇p)*material.k) * dΩ
        end
    end
end

function get_grid()
    # Import grid from abaqus mesh
    grid = get_ferrite_grid(joinpath(@__DIR__, "porous_media", "porous_media_0p75.inp"))

    # Create cellsets for each fieldhandler
    addcellset!(grid, "solid3", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS3")))
    addcellset!(grid, "solid4", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS4R")))
    addcellset!(grid, "porous3", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS3")))
    addcellset!(grid, "porous4", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS4R")))

    # Create faceset for the sides and top
    addfaceset!(grid, "sides", x->(first(x) < eps() || first(x) ≈ 5.0))
    addfaceset!(grid, "top", x->(last(x) ≈ 10.0))
    return grid
end

function create_definition(;t_rise=0.1, p_max=100.0)

    grid = get_grid()

    # Setup the interpolation and integration rules
    dim=Ferrite.getdim(grid)
    ip3_lin = Lagrange{dim, RefTetrahedron, 1}()
    ip4_lin = Lagrange{dim, RefCube, 1}()
    ip3_quad = Lagrange{dim, RefTetrahedron, 2}()
    ip4_quad = Lagrange{dim, RefCube, 2}()
    qr3 = QuadratureRule{dim, RefTetrahedron}(1)
    qr4 = QuadratureRule{dim, RefCube}(2)

    # Setup the MixedDofHandler
    dh = MixedDofHandler(grid)
    push!(dh, FieldHandler([Field(:u, ip3_lin, dim)], getcellset(grid,"solid3")))
    push!(dh, FieldHandler([Field(:u, ip4_lin, dim)], getcellset(grid,"solid4")))
    push!(dh, FieldHandler([Field(:u, ip3_quad, dim), Field(:p, ip3_lin, 1)], getcellset(grid,"porous3")))
    push!(dh, FieldHandler([Field(:u, ip4_quad, dim), Field(:p, ip4_lin, 1)], getcellset(grid,"porous4")))
    close!(dh)

    # Setup cellvalues with the same order as the FieldHandlers in the dh
    # - Linear displacement elements in the solid domain
    # - Taylor hood (quadratic displacement, linear pressure) and linear geometry in porous domain
    cv = ( CellVectorValues(qr3, ip3_lin),
           CellVectorValues(qr4, ip4_lin),
           (u=CellVectorValues(qr3, ip3_quad, ip3_lin), p=CellScalarValues(qr3, ip3_lin)),
           (u=CellVectorValues(qr4, ip4_quad, ip4_lin), p=CellScalarValues(qr4, ip4_lin)) )

    # Add boundary conditions
    # Use `Ferrite.jl` PR427 (temporarily included in FerriteProblems.jl)
    # to make Dirichlet conditions easier and more general
    ch = ConstraintHandler(dh);
    # Fix bottom in y and sides in x
    add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), (x, t) -> zero(Vec{1}), [2]))
    add!(ch, Dirichlet(:u, getfaceset(grid, "sides"), (x,t) -> zero(Vec{1}), [1]))
    # Zero pressure on top surface
    add!(ch, Dirichlet(:p, getfaceset(grid, "top"), (x,t) -> 0.0))
    close!(ch)

    # Add Neumann boundary conditions - normal traction on top
    nh = NeumannHandler(dh);
    add!(nh, Neumann(:u, 2, getfaceset(grid, "top"), (x,t,n) -> -n*clamp(t/t_rise,0,1)*p_max))

    # We then need one material per fieldhandler:
    materials = (Elastic(), Elastic(), PoroElastic(), PoroElastic())

    return FEDefinition(;dh=dh, ch=ch, nh=nh, cv=cv, m=materials, cc=FP.RelativeResidualElementScaling())
end

struct PostProcess{PVD}
    pvd::PVD
    filestem::String
end
function PostProcess(filestem="porous_media")
    pvd = paraview_collection("$filestem.pvd")
    return PostProcess(pvd, filestem)
end

function FESolvers.postprocess!(post::PostProcess, p, step, solver)
    vtk_grid("$(post.filestem)-$step", FP.getdh(p)) do vtk
        vtk_point_data(vtk, FP.getdh(p), FP.getunknowns(p))
        vtk_save(vtk)
        post.pvd[step] = vtk
    end
end

FP.close_postprocessing(post::PostProcess, args...) = vtk_save(post.pvd)

using Logging
Logging.disable_logging(Logging.Warn)
problem = FerriteProblem(create_definition(), PostProcess())
Logging.disable_logging(Logging.LogLevel(-1))
solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper(map(x->x^2, range(0, 1, 41))))
solve_problem!(problem, solver)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

