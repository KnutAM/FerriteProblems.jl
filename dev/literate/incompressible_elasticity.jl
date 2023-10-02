# # Incompressible Elasticity
# This example is the adoption of 
# [`Ferrite.jl`s example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/incompressible_elasticity/)
# To adapt, we make a few changes, specifically 
# that we don't use BlockArrays
# 
#md # The full code, without comments, can be found in the next
#md # [section](@ref incompressible_elasticity-plain-program).

using Ferrite
using FerriteAssembly, FerriteProblems
using FESolvers
import FerriteProblems as FP
import FerriteAssembly as FA

# ## Problem setup
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
function create_values(ipu, ipp)
    RefShape = Ferrite.getrefshape(ipu)
    ## quadrature rules
    qr      = QuadratureRule{RefShape}(3)
    face_qr = FaceQuadratureRule{RefShape}(3)

    ## geometric interpolation
    ip_geo = Lagrange{RefShape,1}()

    ## cell and facevalues for u
    cellvalues_u = CellValues(qr, ipu, ip_geo)
    facevalues_u = FaceValues(face_qr, ipu, ip_geo)

    ## cellvalues for p
    cellvalues_p = CellValues(qr, ipp, ip_geo)

    return cellvalues_u, cellvalues_p, facevalues_u
end;


# We create a DofHandler, with two fields, `:u` and `:p`,
# with possibly different interpolations
function create_dofhandler(grid, ipu, ipp)
    dh = DofHandler(grid)
    add!(dh, :u, ipu) # displacement
    add!(dh, :p, ipp) # pressure
    close!(dh)
    return dh
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

# Define a cache for the material, as used in the original example
function FerriteAssembly.allocate_cell_cache(::LinearElasticity, cv::NamedTuple)
    cellvalues_u = cv[:u]
    return collect([symmetric(shape_gradient(cellvalues_u, 1, i)) for i in 1:getnbasefunctions(cellvalues_u)])
end

function create_definition(ν, ip_u, ip_p)
    grid = create_cook_grid(50, 50)
    dh = create_dofhandler(grid, ip_u, ip_p)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfaceset(dh.grid, "clamped"), Returns(zero(Vec{2}))))
    close!(ch)
    update!(ch, 0.0)

    cv_u, cv_p, fv = create_values(ip_u, ip_p)
    cv = (u=cv_u, p=cv_p)   # Create NamedTuple

    lh = LoadHandler(dh)
    add!(lh, Neumann(:u, fv, getfaceset(grid, "traction"), (x,t,n)->Vec{2}((0.0, 1/16))))

    m = LinearElasticity(;Emod=1.0, ν=ν)
    
    ## Create and return the `FEDefinition`
    domainspec = DomainSpec(dh, m, cv)
    return FEDefinition(domainspec; ch, lh)
end;

# ## Physics (element routine)
# We define the element according to `FerriteAssembly`, and without any Neumann contributions.
# We restrict this element to only work with `FESolvers.LinearProblemSolver` by not including 
# the internal forces in residual. 
# Since the problem results in a symmetric matrix we choose to only assemble the lower part,
# and then symmetrize it after the loop over the quadrature points.
function FerriteAssembly.element_routine!(
    Ke, re, state, ue, mp::LinearElasticity, cv::NamedTuple, buffer
    )
    cellvalues_u = cv[:u]
    cellvalues_p = cv[:p]

    ## Get the local indices for each field
    udofs = dof_range(buffer, :u)
    pdofs = dof_range(buffer, :p)

    ## Extract cached gradients
    ∇Nu_sym_dev = FA.get_user_cache(buffer)

    ## We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    for q_point in 1:getnquadpoints(cellvalues_u)
        for i in 1:length(udofs)
            ∇Nu_sym_dev[i] = dev(symmetric(shape_gradient(cellvalues_u, q_point, i)))
        end
        dΩ = getdetJdV(cellvalues_u, q_point)
        for (i_u, i) in enumerate(udofs)
            for (j_u, j) in enumerate(udofs[1:i_u])
                Ke[i,j] += 2 * mp.G * ∇Nu_sym_dev[i_u] ⊡ ∇Nu_sym_dev[j_u] * dΩ
            end
        end

        for (i_p, i) in enumerate(pdofs)
            δNp = shape_value(cellvalues_p, q_point, i_p)
            for (j_u, j) in enumerate(udofs)
                divδNu = shape_divergence(cellvalues_u, q_point, j_u)
                Ke[i,j] += -δNp * divδNu * dΩ
            end
            for (j_p, j) in enumerate(pdofs[1:i_p])
                Np = shape_value(cellvalues_p, q_point, j_p)
                Ke[i,j] += - 1/mp.K * δNp * Np * dΩ
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
struct IE_PostProcessing
    vtk_file::String
end

function FESolvers.postprocess!(post::IE_PostProcessing, p, step, solver)
    step == 1 && return nothing # We don't want to save the initial conditions. 
    dh = FP.get_dofhandler(p)
    vtk_grid(post.vtk_file, dh) do vtkfile
        vtk_point_data(vtkfile, dh, FP.getunknowns(p))
    end
end;

# ## Solving the problem
function build_problem(;ν, ip_u, ip_p)
    def = create_definition(ν, ip_u, ip_p)
    ip_u_string = isa(ip_u, Lagrange{2,RefTetrahedron,1}) ? "linear" : "quadratic"
    post = IE_PostProcessing("cook_$(ip_u_string)_linear")
    return FerriteProblem(def, post)
end

solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper([0.0,1.0]))
ν = 0.4999999
linear    = Lagrange{2,RefTetrahedron,1}()
quadratic = Lagrange{2,RefTetrahedron,2}()
p1 = build_problem(;ν, ip_u=linear^2,    ip_p=linear)
p2 = build_problem(;ν, ip_u=quadratic^2, ip_p=linear)
solve_problem!(p1, solver)
solve_problem!(p2, solver)

## delete the output                        #src
# Comment out if needed                     #src
rm("cook_linear_linear.vtu"; force=true)    #src
rm("cook_quadratic_linear.vtu"; force=true) #src

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