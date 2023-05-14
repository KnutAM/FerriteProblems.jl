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

function FerriteAssembly.element_residual!(re, state, ae, material::PoroElastic, cv::NamedTuple, buffer)
    # Setup cellvalues and give easier names
    cv_u = cv[:u]
    cv_p = cv[:p]
    num_u = getnbasefunctions(cv_u)
    num_p = getnbasefunctions(cv_p)
    Δt = FerriteAssembly.get_time_increment(buffer)

    # Assign views to the matrix and vector parts
    ae_old = FerriteAssembly.get_aeold(buffer)
    udofs = dof_range(buffer, :u)
    pdofs = dof_range(buffer, :p)
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
end;

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
    fh1 = FieldHandler([Field(:u, ip3_lin, dim)], getcellset(grid,"solid3"))
    add!(dh, fh1)
    fh2 = FieldHandler([Field(:u, ip4_lin, dim)], getcellset(grid,"solid4"))
    add!(dh, fh2)
    fh3 = FieldHandler([Field(:u, ip3_quad, dim), Field(:p, ip3_lin, 1)], getcellset(grid,"porous3"))
    add!(dh, fh3)
    fh4 = FieldHandler([Field(:u, ip4_quad, dim), Field(:p, ip4_lin, 1)], getcellset(grid,"porous4"))
    add!(dh, fh4)
    close!(dh)

    # Setup the AssemblyDomains
    # Solid domain with Triangle elements, linear displacement interpolation
    sdh1 = FerriteAssembly.SubDofHandler(dh, fh1)
    cv1 = CellVectorValues(qr3, ip3_lin)
    ad1 = AssemblyDomain("solid3", sdh1, Elastic(), cv1)

    # Solid domain with Quadrilateral elements, linear displacement interpolation
    sdh2 = FerriteAssembly.SubDofHandler(dh, fh2)
    cv2 = CellVectorValues(qr4, ip4_lin)
    ad2 = AssemblyDomain("solid4", sdh2, Elastic(), cv2)

    # Porous domain with Triangle elements
    # Taylor hood: (quadratic displacement and linear pressure interpolation)
    sdh3 = FerriteAssembly.SubDofHandler(dh, fh3)
    cv3 = (u=CellVectorValues(qr3, ip3_quad, ip3_lin), p=CellScalarValues(qr3, ip3_lin))
    ad3 = AssemblyDomain("porous3", sdh3, PoroElastic(), cv3)

    # Porous domain with Quadrilateral elements
    # Taylor hood: (quadratic displacement and linear pressure interpolation)
    sdh4 = FerriteAssembly.SubDofHandler(dh, fh4)
    cv4 = (u=CellVectorValues(qr4, ip4_quad, ip4_lin), p=CellScalarValues(qr4, ip4_lin))
    ad4 = AssemblyDomain("porous4", sdh4, PoroElastic(), cv4)

    # Add boundary conditions
    ch = ConstraintHandler(dh);
    # Fix bottom in y and sides in x
    add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), Returns(0.0), [2]))
    add!(ch, Dirichlet(:u, getfaceset(grid, "sides"), Returns(0.0), [1]))
    # Zero pressure on top surface
    add!(ch, Dirichlet(:p, getfaceset(grid, "top"), Returns(0.0)))
    close!(ch)

    # Add Neumann boundary conditions - normal traction on top
    nh = NeumannHandler(dh);
    add!(nh, Neumann(:u, 2, getfaceset(grid, "top"), (x,t,n) -> -n*clamp(t/t_rise,0,1)*p_max))

    return FEDefinition([ad1, ad2, ad3, ad4]; ch=ch, nh=nh)
end;

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

FP.close_postprocessing(post::PostProcess, args...) = vtk_save(post.pvd);

problem = FerriteProblem(create_definition(), PostProcess())
solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper(map(x->x^2, range(0, 1, 41))))
solve_problem!(problem, solver)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

