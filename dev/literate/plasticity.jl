# # Plasticity

# This example is taken from 
# [Ferrite.jl's plasticity example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/)
# and shows how `FerriteProblems` can be used to simplify the setup of this nonlinear problem with 
# time dependent loading.
# 
# First we need to load all required packages
using Ferrite, Tensors, SparseArrays, LinearAlgebra
using FerriteProblems, FESolvers, FerriteNeumann, FerriteAssembly
import FerriteProblems as FP
using Plots; gr()

# We then define the material by including the definitions used in the original 
# [example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/),
# by using the [J2Plasticity.jl file](J2Plasticity.jl). We have modified the 
# names and functions to comply with 
# [`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)
include("J2Plasticity.jl");

# This file defines the 
# * `J2Plasticity<:AbstractMaterial` material type with the constructor `J2Plasticity(E,ν,σy0,H)`
# * State variable struct `J2PlasticityMaterialState<:AbstractMaterialState`
# * The MaterialModelsBase function `initial_material_state` and `material_response`
# For a single material according to the `MaterialModelsBase` interface, the element is already 
# included in FerriteProblems, see `src/MaterialModelsBase.jl`

# ## Problem definition
# We first create the problem's definition

traction_function(time) = time*1.e7 # N/m² 

function setup_problem_definition()
    ## Define material properties ("J2Plasticity.jl" file)
    material = J2Plasticity(200.0e9, 0.3, 200.e6, 10.0e9)
    
    ## Cell and facevalues (`Ferrite.jl`)
    interpolation = Lagrange{3, RefTetrahedron, 1}()
    cv = CellVectorValues(QuadratureRule{3,RefTetrahedron}(2), interpolation)
    fv = FaceVectorValues(QuadratureRule{2,RefTetrahedron}(3), interpolation)

    ## Grid and degrees of freedom (`Ferrite.jl`)
    grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.,1.,1.)))
    dh = DofHandler(grid); push!(dh, :u, 3, interpolation); close!(dh)

    ## Constraints (Dirichlet boundary conditions, `Ferrite.jl`)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfaceset(grid, "left"), (x,t) -> zeros(3), [1, 2, 3]))
    close!(ch)

    ## Neumann boundary conditions (`FerriteNeumann.jl`)
    nh = NeumannHandler(dh)
    add!(nh, Neumann(:u, fv, getfaceset(grid, "right"), (x,t,n)->Vec{3}((0.0, 0.0, traction_function(t)))))

    ## Initial material states (using FerriteAssembly's `create_states`)
    states = create_states(dh, x->initial_material_state(material), cv)

    return FEDefinition(;dh=dh, ch=ch, nh=nh, cv=cv, m=material, initialstate=states)
end;

# For the problem at hand, we need to define the element routine, 
# following `FerriteAssembly`s interface. This function is almost equivalent to 
# the `assemble_cell!` in the original example, except that 
# 1) We don't have to `reinit!` as `FerriteAssembly` does that before calling
# 2) The traction is handled by `FerriteNeumann` and is not done for each cell
# For convenience, this element is already implemented for materials of 
# MaterialModelsBase type `AbstractMaterial`, and we don't need to define it here. 

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
# Note that, internally, FerriteProblems imports the FESolvers functions 
# `getunknowns`, `getjacobian`, and `getresidual`, such that you can access these 
# via `FerriteProblems.` (or `FP.` if using the `import FerriteProblems as FP` above).
# For convenience, `FerriteProblems` will call `FESolvers.postprocess!` with the 
# `post` as the first argument making it easy to dispatch on: 
function FESolvers.postprocess!(post::PlasticityPostProcess, p, step, solver)
    ## p::FerriteProblem
    ## First, we save some values directly in the `post` struct
    push!(post.tmag, traction_function(FP.gettime(p)))
    push!(post.umag, maximum(abs, FP.getunknowns(p)))

    ## Second, we save some results to file
    ## * We must always start by adding the next step.
    FP.addstep!(p.io, p)
    ## * Save the dof values (only displacments in this case)
    FP.savedofdata!(p.io, FP.getunknowns(p))
    ## * Save the state in each integration point
    FP.saveipdata!(p.io, FP.getstate(p), "state")
end;

# We also define a helper function to plot the results after completion
function plot_results(problem::FerriteProblem{<:PlasticityPostProcess}; 
    plt=plot(), label=nothing, markershape=:auto, markersize=4
    )
    umax = vcat(0.0, problem.post.umag)
    tmag = vcat(0.0, problem.post.tmag)
    plot!(plt, umax, tmag, linewidth=0.5, title="Traction-displacement", label=label, 
        markeralpha=0.75, markershape=markershape, markersize=markersize)
    ylabel!(plt, "Traction [Pa]")
    xlabel!(plt, "Maximum deflection [m]")
    return plt
end;

# ## Solving the problem
# Finally, we can solve the problem with different time stepping strategies 
# and plot the results. Here, we use `FerriteProblems`' `safesolve` that 
# (1) creates our full `problem::FerriteProblem` 
# and (2) ensures that files are closed even when the problem doesn't converge. 

global umax_solution = [0.0] # To save result for test #hide

function example_solution()
    def = setup_problem_definition()

    ## Fixed uniform time steps
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(;num_steps=25,Δt=0.04))
    problem = FerriteProblem(def, PlasticityPostProcess(), joinpath(pwd(), "A"))
    solve_problem!(solver, problem)
    plt = plot_results(problem, label="uniform", markershape=:x, markersize=5)

    ## Same time steps as Ferrite example, overwrite results by specifying the same folder
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(append!([0.], collect(0.5:0.05:1.0))))
    problem = FerriteProblem(def, PlasticityPostProcess(), joinpath(pwd(), "A"))
    solve_problem!(solver, problem)
    plot_results(problem, plt=plt, label="fixed", markershape=:circle)
    umax_solution[1] = problem.post.umag[end] # Save value for comparison  #hide

    ## Adaptive time stepping, save results to new folder
    ts = AdaptiveTimeStepper(0.05, 1.0; Δt_min=0.01, Δt_max=0.2)
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0, maxiter=6), ts)
    problem = FerriteProblem(def, PlasticityPostProcess(), joinpath(pwd(), "B"))
    solve_problem!(solver, problem)
    plot_results(problem, plt=plt, label="adaptive", markershape=:circle)
    
    plot!(;legend=:bottomright)
    return plt, problem, solver
end;

plt, problem, solver = example_solution();

using Test # Compare to Ferrite.jl's example #hide
@test isapprox(umax_solution[1], 0.254452; rtol=1.e-4);  #hide

# Which gives the following result when running `display(plt)`
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