# # Plasticity

# This example is taken from 
# [Ferrite.jl's plasticity example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/)
# and shows how `FerriteProblems` can be used to simplify the setup of this nonlinear problem with 
# time dependent loading.
# 
# First we need to load all required packages
using Ferrite, Tensors, SparseArrays, LinearAlgebra, Plots
using FerriteProblems, FESolvers, FerriteAssembly, FerriteNeumann

# We then define the material by including the definitions used in the original 
# [example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/),
# by using the [J2Plasticity.jl file](J2Plasticity.jl)
include("J2Plasticity.jl");

# ## Problem definition
# We first create the problem's definition

traction_function(time) = time*1.e7 # N/m²

function setup_problem_definition()
    ## Define material properties
    E = 200.0e9
    material = J2Plasticity(E, 0.3, 200.e6, E/20)
    
    ## Mesh
    L = 10.0; w = 1.0; h = 1.0  # Dimensions
    n = 2
    nels = (10n, n, 2n)
    P1 = Vec((0.0, 0.0, 0.0))
    P2 = Vec((L, w, h))
    grid = generate_grid(Tetrahedron, nels, P1, P2)
    
    ## Interpolation and integration => FEValues
    interpolation = Lagrange{3, RefTetrahedron, 1}()
    ## setup quadrature rules
    qr      = QuadratureRule{3,RefTetrahedron}(2)
    face_qr = QuadratureRule{2,RefTetrahedron}(3)

    ## create geometric interpolation (use the same as for displacements)
    interpolation_geom = Lagrange{3,RefTetrahedron,1}()

    ## cell and facevalues
    cv = CellVectorValues(qr, interpolation, interpolation_geom)
    fv = FaceVectorValues(face_qr, interpolation, interpolation_geom)

    ## Degrees of freedom
    dh = DofHandler(grid)
    push!(dh, :u, 3, interpolation) # add a displacement field with 3 components
    close!(dh)

    ## Constraints (Dirichlet boundary conditions)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfaceset(grid, "left"), (x,t) -> [0.0, 0.0, 0.0], [1, 2, 3]))
    close!(ch)

    ## Neumann boundary conditions
    nh = NeumannHandler(dh)
    add!(nh, Neumann(:u, fv, getfaceset(grid, "right"), (x,t,n)->Vec{3}((0.0, 0.0, traction_function(t)))))

    ## Initial material states
    states = [ [J2PlasticityMaterialState() for _ in 1:getnquadpoints(cv)] for _ in 1:getncells(grid)]
    

    return FEDefinition(dh, ch, nh, cv, material, states)
end;

# For the problem at hand, we need to define the element routine, 
# following `FerriteAssembly`s interface. This function is almost equivalent to 
# the `assemble_cell!` in the original example, except that 
# 1) We don't have to `reinit!` as `FerriteAssembly` does that before calling
# 2) The traction is handled by `FerriteNeumann` and is not done for each cell

function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, 
    ue::AbstractVector, ae_old::AbstractVector,
    state::AbstractVector, material::J2Plasticity, cellvalues::CellVectorValues, 
    dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
    )
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        ## For each integration point, compute stress and material stiffness
        ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
        σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state[q_point])

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
            re[i] += (δϵ ⊡ σ) * dΩ # add internal force to residual
            for j in 1:i # loop only over lower half
                Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke)
end;

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

# At this point, we can define 
# `problem = FerriteProblem(setup_problem_definition())`, which can be solved 
# with `FESolvers.jl`'s `solve_problem!`. 
# But to get any results apart from the final state, we need to define 
# the postprocessing after each step.

# ## Setup postprocessing
# In contrast to the original example, 
# we do not save directly to a vtk-file, but use `FerriteProblem`'s IO features 
# to save to JLD2 files. This has the advantage that further postprocessing can 
# be done after the simulation, and we can then choose to export to the VTK-format 
# or plot directly using e.g. `FerriteViz.jl`. 
# We start by defining our custom postprocessing type.
struct PlasticityPostProcess{T}
    tmag::Vector{T}
    umag::Vector{T}
end
PlasticityPostProcess() = PlasticityPostProcess(Float64[], Float64[]);

# With this postprocessing type, we can now define the postprocessing in FESolvers.
# Since the `FerriteProblem` type is used for all types of problems, it can be 
# useful to dispatch on the contained postprocessing type to have multiple problems 
# defined in the same session/package. 
function FESolvers.postprocess!(p::FerriteProblem{<:PlasticityPostProcess}, step, solver)
    ## First, we save some values directly in the `post` struct 
    push!(p.post.tmag, traction_function(FerriteProblems.gettime(p)))
    push!(p.post.umag, maximum(abs, FESolvers.getunknowns(p)))

    ## Second, we save some results to file
    ## * We must always start by adding the next step. 
    FerriteProblems.addstep!(p.io, FerriteProblems.gettime(p))
    ## * Save the dof values (only displacments in this case)
    FerriteProblems.savedofdata!(p.io, FESolvers.getunknowns(p))
    ## * Save the state in each integration point
    FerriteProblems.saveipdata!(p.io, FerriteProblems.getstate(p), "state")
end;

# ## Solving the problem
# First, we define a helper function to plot the results after the solution
function plot_results(problem::FerriteProblem{<:PlasticityPostProcess}; plt=plot(), label=nothing, markershape=:auto, markersize=4)
    umax = vcat(0.0, problem.post.umag)
    tmag = vcat(0.0, problem.post.tmag)
    plot!(plt, umax, tmag, linewidth=0.5, title="Traction-displacement", label=label, 
        markeralpha=0.75, markershape=markershape, markersize=markersize)
    ylabel!(plt, "Traction [Pa]")
    xlabel!(plt, "Maximum deflection [m]")
    return plt
end;

# **Temporary:** 
# Wrap `solve_problem!` to make it save any open stuff at the end, 
# even in case of an exception.
# This is temporary, but would be nice to handle within FESolvers
function wrapped_solve!(solver, problem)
    try
        solve_problem!(solver, problem)
    finally
        FerriteProblems.close_problem(problem)
    end
end;

# Finally, we can solve the problem with different time stepping strategies 
# and plot the results
function example_solution()
    def = setup_problem_definition()
    makeproblem(_def, folder) = FerriteProblem(_def, PlasticityPostProcess(), joinpath(pwd(), folder))

    ## Fixed uniform time steps
    problem = makeproblem(def, "A")
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(;num_steps=25,Δt=0.04))
    wrapped_solve!(solver, problem)
        
    plt = plot_results(problem, label="uniform", markershape=:x, markersize=5)

    ## Same time steps as Ferrite example
    problem = makeproblem(def, "A")
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(append!([0.], collect(0.5:0.05:1.0))))
    wrapped_solve!(solver, problem)
    plot_results(problem, plt=plt, label="fixed", markershape=:circle)

    ## Adaptive time stepping 
    problem = makeproblem(def, "C")
    ts = AdaptiveTimeStepper(0.05, 1.0; Δt_min=0.01, Δt_max=0.2)
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0, maxiter=4), ts)
    wrapped_solve!(solver, problem)

    plot_results(problem, plt=plt, label="adaptive", markershape=:circle)
    plot!(;legend=:bottomright)
    return plt, problem, solver
end;

plt, problem, solver = example_solution();
display(plt);

# Which gives the following result
# 
# 
# ![](plasticity.svg)

#md # ## Plain program
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`plasticity.jl`](plasticity.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```