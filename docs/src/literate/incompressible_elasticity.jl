# # Incompressible Elasticity
# This example is the adoption of 
# [`Ferrite.jl`s example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/incompressible_elasticity/)
# To adapt, we make a few changes, specifically 
# that we don't use BlockArrays
# 
#md # The full code, without comments, can be found in the next
#md # [section](@ref incompressible_elasticity-plain-program).

using Ferrite
using FerriteNeumann, FerriteAssembly, FerriteProblems
using FESolvers
import FerriteProblems as FP

## Problem setup
# First we generate a simple grid, specifying the 4 corners of Cooks membrane.
function create_cook_grid(nx, ny)
    corners = [Vec{2}((0.0,   0.0)),
               Vec{2}((48.0, 44.0)),
               Vec{2}((48.0, 60.0)),
               Vec{2}((0.0,  44.0))]
    grid = generate_grid(Triangle, (nx, ny), corners);
    ## facesets for boundary conditions
    addfaceset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0);
    addfaceset!(grid, "traction", x -> norm(x[1]) ≈ 48.0);
    return grid
end;

# Next we define a function to set up our cell- and facevalues.
function create_values(interpolation_u, interpolation_p)
    ## quadrature rules
    qr      = QuadratureRule{2,RefTetrahedron}(3)
    face_qr = QuadratureRule{1,RefTetrahedron}(3)

    ## geometric interpolation
    interpolation_geom = Lagrange{2,RefTetrahedron,1}()

    ## cell and facevalues for u
    cellvalues_u = CellVectorValues(qr, interpolation_u, interpolation_geom)
    facevalues_u = FaceVectorValues(face_qr, interpolation_u, interpolation_geom)

    ## cellvalues for p
    cellvalues_p = CellScalarValues(qr, interpolation_p, interpolation_geom)

    return cellvalues_u, cellvalues_p, facevalues_u
end;


# We create a DofHandler, with two fields, `:u` and `:p`,
# with possibly different interpolations
function create_dofhandler(grid, ipu, ipp)
    dh = DofHandler(grid)
    push!(dh, :u, 2, ipu) # displacement
    push!(dh, :p, 1, ipp) # pressure
    close!(dh)
    return dh
end;

# We also need to add Dirichlet boundary conditions on the `"clamped"` faceset.
# We specify a homogeneous Dirichlet bc on the displacement field, `:u`.
function create_bc(dh)
    dbc = ConstraintHandler(dh)
    add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "clamped"), (x,t) -> zero(Vec{2}), [1,2]))
    close!(dbc)
    t = 0.0
    update!(dbc, t)
    return dbc
end;

# The material is linear elastic, which is here specified by the shear and bulk moduli
struct LinearElasticity{T}
    G::T
    K::T
end
function LinearElasticity(;Emod, ν)
    Gmod = Emod / 2(1 + ν)
    Kmod = Emod * ν / ((1+ν) * (1-2ν))
    return LinearElasticity(Gmod, Kmod)
end

function create_definition(ν, ip_u, ip_p)
    grid = create_cook_grid(50, 50)
    dh = create_dofhandler(grid, ip_u, ip_p)
    ch = create_bc(dh)
    cv_u, cv_p, fv = create_values(ip_u, ip_p)
    cv = (u=cv_u, p=cv_p)   # Create NamedTuple

    nh = NeumannHandler(dh)
    add!(nh, Neumann(:u, fv, getfaceset(grid, "traction"), (x,t,n)->Vec{2}((0.0, 1/16))))

    m = LinearElasticity(;Emod=1.0, ν=ν)
    
    ## Create and return the `FEDefinition`
    return FEDefinition(;dh=dh, ch=ch, nh=nh, cv=cv, m=m)
end;

# ## Physics (element routine)
# We define the element according to `FerriteAssembly`, and without any Neumann contributions.
# We restrict this element to only work with `FESolvers.LinearProblemSolver` by not including 
# the internal forces in residual. 
# Since the problem results in a symmetric matrix we choose to only assemble the lower part,
# and then symmetrize it after the loop over the quadrature points.
function FerriteAssembly.element_routine!(
    Ke, re, state, ue, mp::LinearElasticity, cv::NamedTuple, dh_fh, Δt, buffer
    )
    cellvalues_u = cv[:u]
    cellvalues_p = cv[:p]
    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    ## Create views of the different fields
    udofs = dof_range(dh_fh, :u)
    pdofs = dof_range(dh_fh, :p)
    Ke_uu = @view Ke[udofs,udofs]
    Ke_pu = @view Ke[pdofs,udofs]
    Ke_pp = @view Ke[pdofs,pdofs]

    ## Temporary, should be cached
    ∇Nu_sym_dev = zeros(SymmetricTensor{2,2,Float64,3}, n_basefuncs_u)

    ## We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    for q_point in 1:getnquadpoints(cellvalues_u)
        for i in 1:n_basefuncs_u
            ∇Nu_sym_dev[i] = dev(symmetric(shape_gradient(cellvalues_u, q_point, i)))
        end
        dΩ = getdetJdV(cellvalues_u, q_point)
        for i in 1:n_basefuncs_u
            for j in 1:i
                Ke_uu[i,j] += 2 * mp.G * ∇Nu_sym_dev[i] ⊡ ∇Nu_sym_dev[j] * dΩ
            end
        end

        for i in 1:n_basefuncs_p
            δNp = shape_value(cellvalues_p, q_point, i)
            for j in 1:n_basefuncs_u
                divδNu = shape_divergence(cellvalues_u, q_point, j)
                Ke_pu[i,j] += -δNp * divδNu * dΩ
            end
            for j in 1:i
                Np = shape_value(cellvalues_p, q_point, j)
                Ke_pp[i,j] += - 1/mp.K * δNp * Np * dΩ
            end
        end
    end
    symmetrize_lower!(Ke)
end

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

# ## Postprocessing
struct PostProcessing
    vtk_file::String
end

function FESolvers.postprocess!(post::PostProcessing, p, step, solver)
    dh = FP.getdh(p)
    vtk_grid(post.vtk_file, dh) do vtkfile
        vtk_point_data(vtkfile, dh, FP.getunknowns(p))
    end
end;

# ## Solving the problem
function build_problem(ν, ip_u, ip_p)
    def = create_definition(ν, ip_u, ip_p)
    ip_u_string = isa(ip_u, Lagrange{2,RefTetrahedron,1}) ? "linear" : "quadratic"
    post = PostProcessing("cook_$(ip_u_string)_linear")
    return FerriteProblem(def, post)
end

solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper([0.0,1.0]))
ν = 0.4999999
linear    = Lagrange{2,RefTetrahedron,1}()
quadratic = Lagrange{2,RefTetrahedron,2}()
p1 = build_problem(ν, linear, linear)
p2 = build_problem(ν, quadratic, linear)
solve_problem!(p1, solver)
solve_problem!(p2, solver)

## test the result                 #src
using Test                         #src
@test norm(FP.getunknowns(p2)) ≈ 919.2122668839389 #src

#md # ## [Plain program](@id incompressible_elasticity-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here:
#md # [`incompressible_elasticity.jl`](incompressible_elasticity.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```