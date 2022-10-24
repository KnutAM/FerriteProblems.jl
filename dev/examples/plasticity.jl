using Ferrite, Tensors, SparseArrays, LinearAlgebra
using FerriteProblems, FESolvers, FerriteAssembly, FerriteNeumann

using Plots; gr()

include("J2Plasticity.jl");

traction_function(time) = time*1.e7 # N/m²

function setup_problem_definition()
    # Define material properties
    material = J2Plasticity(200.0e9, 0.3, 200.e6, 10.0e9)

    # Cell and facevalues
    interpolation = Lagrange{3, RefTetrahedron, 1}()
    cv = CellVectorValues(QuadratureRule{3,RefTetrahedron}(2), interpolation)
    fv = FaceVectorValues(QuadratureRule{2,RefTetrahedron}(3), interpolation)

    # Grid and degrees of freedom
    grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.,1.,1.)))
    dh = DofHandler(grid); push!(dh, :u, 3, interpolation); close!(dh)

    # Constraints (Dirichlet boundary conditions)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfaceset(grid, "left"), (x,t) -> zeros(3), [1, 2, 3]))
    close!(ch)

    # Neumann boundary conditions
    nh = NeumannHandler(dh)
    add!(nh, Neumann(:u, fv, getfaceset(grid, "right"), (x,t,n)->Vec{3}((0.0, 0.0, traction_function(t)))))

    # Initial material states
    states = create_states(dh, x->J2PlasticityMaterialState(), cv)

    return FEDefinition(;dh=dh, ch=ch, nh=nh, cv=cv, m=material, initialstate=states)
end;

function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, state::AbstractVector,
    ue::AbstractVector, material::J2Plasticity, cellvalues::CellVectorValues,
    args...)
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        # For each integration point, compute stress and material stiffness
        ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
        σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state[q_point])

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
            re[i] += (δϵ ⊡ σ) * dΩ # add internal force to residual
            for j in 1:n_basefuncs
                Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
            end
        end
    end
end;

struct PlasticityPostProcess{T}
    tmag::Vector{T}
    umag::Vector{T}
end
PlasticityPostProcess() = PlasticityPostProcess(Float64[], Float64[]);

function FESolvers.postprocess!(post::PlasticityPostProcess, p, step, solver)
    # p::FerriteProblem
    # First, we save some values directly in the `post` struct
    push!(post.tmag, traction_function(FerriteProblems.gettime(p)))
    push!(post.umag, maximum(abs, FESolvers.getunknowns(p)))

    # Second, we save some results to file
    # * We must always start by adding the next step.
    FerriteProblems.addstep!(p.io, FerriteProblems.gettime(p))
    # * Save the dof values (only displacments in this case)
    FerriteProblems.savedofdata!(p.io, FESolvers.getunknowns(p))
    # * Save the state in each integration point
    FerriteProblems.saveipdata!(p.io, FerriteProblems.getstate(p), "state")
end;

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

global umax_solution = [0.0] # To save result for test #hide

function example_solution()
    def = setup_problem_definition()

    # Fixed uniform time steps
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(;num_steps=25,Δt=0.04))
    problem = safesolve(solver, def, PlasticityPostProcess(), joinpath(pwd(), "A"))
    plt = plot_results(problem, label="uniform", markershape=:x, markersize=5)

    # Same time steps as Ferrite example, overwrite results by specifying the same folder
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(append!([0.], collect(0.5:0.05:1.0))))
    problem = safesolve(solver, def, PlasticityPostProcess(), joinpath(pwd(), "A"))
    plot_results(problem, plt=plt, label="fixed", markershape=:circle)
    umax_solution[1] = problem.post.umag[end] # Save value for comparison  #hide

    # Adaptive time stepping, save results to new folder
    ts = AdaptiveTimeStepper(0.05, 1.0; Δt_min=0.01, Δt_max=0.2)
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0, maxiter=6), ts)
    problem = safesolve(solver, def, PlasticityPostProcess(), joinpath(pwd(), "B"))
    plot_results(problem, plt=plt, label="adaptive", markershape=:circle)

    plot!(;legend=:bottomright)
    return plt, problem, solver
end;

plt, problem, solver = example_solution();

using Test # Compare to Ferrite.jl's example #hide
@test isapprox(umax_solution[1], 0.254452; rtol=1.e-4);  #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

