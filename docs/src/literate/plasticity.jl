# # Plasticity

# This example is taken from 
# [Ferrite.jl's plasticity example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/)
# and shows how `FerriteProblems` can be used to simplify the setup of this nonlinear problem with 
# time dependent loading.
# 
# First we need to load all required packages
using Ferrite, Tensors, SparseArrays, LinearAlgebra
using FerriteProblems, FESolvers, FerriteAssembly
import FerriteProblems as FP
import FerriteAssembly.ExampleElements: J2Plasticity
using Plots; gr();

# The material from the original 
# [example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/),
# is available in `FerriteAssembly.ExampleElements` as `J2Plasticity`, 
# serving as an example of a material that complies with the  
# [`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)
# interface.

# ## Problem definition
# We first create the problem's definition. 
# To be able to save the results using JLD2, we cannot use anonymous functions,
# so a struct for the ramping is created instead:
struct VectorRamp{dim,T}<:Function
    ramp::Vec{dim,T}
end
(vr::VectorRamp)(x, t, n) = t*vr.ramp
const traction_function = VectorRamp(Vec(0.0, 0.0, 1.e7))
function setup_problem_definition()
    ## Define material properties
    material = J2Plasticity(;E=200.0e9, ν=0.3, σ0=200.e6, H=10.0e9)
    ip =  Lagrange{RefTetrahedron, 1}()^3
    ## CellValues
    cv = CellValues(QuadratureRule{RefTetrahedron}(2), ip)

    ## Grid and degrees of freedom (`Ferrite.jl`)
    grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.,1.,1.)))
    dh = DofHandler(grid); push!(dh, :u, ip); close!(dh)

    ## Constraints (Dirichlet boundary conditions, `Ferrite.jl`)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfaceset(grid, "left"), Returns(zero(Vec{3}))))
    close!(ch)

    ## Neumann boundary conditions
    lh = LoadHandler(dh)
    quad_order = 3
    add!(lh, Neumann(:u, quad_order, getfaceset(grid, "right"), traction_function))

    domainspec = DomainSpec(dh, material, cv)
    return FEDefinition(domainspec; ch, lh)
end;

# For the problem at hand, `FerriteAssembly.element_routine!` is defined in `FerriteAssembly.jl`.

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
function FESolvers.postprocess!(post::PlasticityPostProcess, p, solver)
    ## p::FerriteProblem
    ## First, we save some values directly in the `post` struct
    push!(post.tmag, traction_function(zero(Vec{3}), FP.get_time(p), zero(Vec{3}))[3])
    push!(post.umag, maximum(abs, FP.getunknowns(p)))

    ## Second, we save some results to file
    ## * We must always start by adding the next step.
    FP.addstep!(p.io, p)
    ## * Save the dof values (only displacments in this case)
    FP.savedofdata!(p.io, FP.getunknowns(p))
    ## * Save the state in each integration point
    FP.saveipdata!(p.io, FP.get_state(p), "state")
end;

# We also define a helper function to plot the results after completion
function plot_results(post::PlasticityPostProcess; 
    plt=plot(), label=nothing, markershape=:auto, markersize=4
    )
    plot!(plt, post.umag, post.tmag, linewidth=0.5, title="Traction-displacement", label=label, 
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
    solve_problem!(problem, solver)
    plt = plot_results(problem.post, label="uniform", markershape=:x, markersize=5)

    ## Same time steps as Ferrite example, overwrite results by specifying the same folder
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(append!([0.], collect(0.5:0.05:1.0))))
    problem = FerriteProblem(def, PlasticityPostProcess(), joinpath(pwd(), "A"))
    solve_problem!(problem, solver)
    plot_results(problem.post, plt=plt, label="fixed", markershape=:circle)
    umax_solution[1] = problem.post.umag[end] # Save value for comparison  #hide

    ## Adaptive time stepping, save results to new folder
    ts = AdaptiveTimeStepper(0.05, 1.0; Δt_min=0.01, Δt_max=0.2)
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0, maxiter=6), ts)
    problem = FerriteProblem(def, PlasticityPostProcess(), joinpath(pwd(), "B"))
    solve_problem!(problem, solver)
    plot_results(problem.post, plt=plt, label="adaptive", markershape=:circle)
    
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