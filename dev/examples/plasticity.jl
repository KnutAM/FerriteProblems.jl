using Ferrite, Tensors, SparseArrays, LinearAlgebra
using FerriteProblems, FESolvers, FerriteAssembly, FerriteNeumann

using Plots; gr()

include("J2Plasticity.jl");

traction_function(time) = time*1.e7 # N/m²

function setup_problem_definition()
    # Define material properties
    material = J2Plasticity(200.0e9, 0.3, 200.e6, 10.0e9)

    # Mesh
    grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.,1.,1.)))

    # Cell and facevalues
    interpolation = Lagrange{3, RefTetrahedron, 1}()
    cv = CellVectorValues(QuadratureRule{3,RefTetrahedron}(2), interpolation)
    fv = FaceVectorValues(QuadratureRule{2,RefTetrahedron}(3), interpolation)

    # Degrees of freedom
    dh = DofHandler(grid)
    push!(dh, :u, 3, interpolation) # add a displacement field with 3 components
    close!(dh)

    # Constraints (Dirichlet boundary conditions)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfaceset(grid, "left"), (x,t) -> [0.0, 0.0, 0.0], [1, 2, 3]))
    close!(ch)

    # Neumann boundary conditions
    nh = NeumannHandler(dh)
    add!(nh, Neumann(:u, fv, getfaceset(grid, "right"), (x,t,n)->Vec{3}((0.0, 0.0, traction_function(t)))))

    # Initial material states
    states = [ [J2PlasticityMaterialState() for _ in 1:getnquadpoints(cv)] for _ in 1:getncells(grid)]

    return FEDefinition(dh, ch, nh, cv, material, states)
end;

function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector,
    ue::AbstractVector, ae_old::AbstractVector,
    state::AbstractVector, material::J2Plasticity, cellvalues::CellVectorValues,
    dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
    )
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        # For each integration point, compute stress and material stiffness
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

struct PlasticityPostProcess{T}
    tmag::Vector{T}
    umag::Vector{T}
end
PlasticityPostProcess() = PlasticityPostProcess(Float64[], Float64[]);

function FESolvers.postprocess!(p::FerriteProblem{<:PlasticityPostProcess}, step, solver)
    # First, we save some values directly in the `post` struct
    push!(p.post.tmag, traction_function(FerriteProblems.gettime(p)))
    push!(p.post.umag, maximum(abs, FESolvers.getunknowns(p)))

    # Second, we save some results to file
    # * We must always start by adding the next step.
    FerriteProblems.addstep!(p.io, FerriteProblems.gettime(p))
    # * Save the dof values (only displacments in this case)
    FerriteProblems.savedofdata!(p.io, FESolvers.getunknowns(p))
    # * Save the state in each integration point
    FerriteProblems.saveipdata!(p.io, FerriteProblems.getstate(p), "state")
end;

function plot_results(problem::FerriteProblem{<:PlasticityPostProcess}; plt=plot(), label=nothing, markershape=:auto, markersize=4)
    umax = vcat(0.0, problem.post.umag)
    tmag = vcat(0.0, problem.post.tmag)
    plot!(plt, umax, tmag, linewidth=0.5, title="Traction-displacement", label=label,
        markeralpha=0.75, markershape=markershape, markersize=markersize)
    ylabel!(plt, "Traction [Pa]")
    xlabel!(plt, "Maximum deflection [m]")
    return plt
end;

function wrapped_solve!(solver, problem)
    try
        solve_problem!(solver, problem)
    finally
        FerriteProblems.close_problem(problem)
    end
end;

function example_solution()
    def = setup_problem_definition()
    makeproblem(_def, folder) = FerriteProblem(_def, PlasticityPostProcess(), joinpath(pwd(), folder))

    # Fixed uniform time steps
    problem = makeproblem(def, "A")
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(;num_steps=25,Δt=0.04))
    wrapped_solve!(solver, problem)

    plt = plot_results(problem, label="uniform", markershape=:x, markersize=5)

    # Same time steps as Ferrite example, overwrite results
    problem = makeproblem(def, "A")
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0), FixedTimeStepper(append!([0.], collect(0.5:0.05:1.0))))
    wrapped_solve!(solver, problem)
    plot_results(problem, plt=plt, label="fixed", markershape=:circle)

    # Adaptive time stepping, save results to new folder
    problem = makeproblem(def, "B")
    ts = AdaptiveTimeStepper(0.05, 1.0; Δt_min=0.01, Δt_max=0.2)
    solver = QuasiStaticSolver(NewtonSolver(;tolerance=1.0, maxiter=6), ts)
    wrapped_solve!(solver, problem)

    plot_results(problem, plt=plt, label="adaptive", markershape=:circle)
    plot!(;legend=:bottomright)
    return plt, problem, solver
end;

plt, problem, solver = example_solution();

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
