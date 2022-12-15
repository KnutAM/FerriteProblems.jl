using Ferrite, FerriteProblems, FerriteAssembly, FESolvers
import FerriteProblems as FP

Base.@kwdef struct FicksLaw{T}
    k::T=1.0e-3    # Thermal conductivity
    f::T=5.0e-1    # Constant heat source
end

function FerriteAssembly.element_routine!(
    Ke, re, state, ue, m::FicksLaw, cellvalues, dh_fh, Δt, buffer
    )
    ue_old = buffer.ae_old  # Extract old values from the CellBuffer (TODO: Perhaps good to have get-functions for this)
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        u = function_value(cellvalues, q_point, ue)
        uold = function_value(cellvalues, q_point, ue_old)
        ∇u = function_gradient(cellvalues, q_point, ue)
        for i in 1:n_basefuncs
            δN = shape_value(cellvalues, q_point, i)
            ∇δN = shape_gradient(cellvalues, q_point, i)
            re[i] += (δN * (u - uold - Δt * m.f) + Δt * m.k * ∇δN ⋅ ∇u) * dΩ
            for j in 1:n_basefuncs
                N = shape_value(cellvalues, q_point, j)
                ∇N = shape_gradient(cellvalues, q_point, j)
                Ke[i, j] += (δN*N + Δt * m.k * (∇δN ⋅ ∇N)) * dΩ
            end
        end
    end
end;

function create_definition()
    # **Grid**
    grid = generate_grid(Quadrilateral, (100, 100));

    # **Cell values**
    cellvalues = CellScalarValues(
        QuadratureRule{2, RefCube}(2),
        Lagrange{2, RefCube, 1}());

    # **Degrees of freedom**
    # After this, we can define the `DofHandler` and distribute the DOFs of the problem.
    dh = DofHandler(grid); push!(dh, :u, 1); close!(dh)

    # **Boundary conditions**
    # Zero pressure on $\partial \Omega_1$ and linear ramp followed by constant pressure on $\partial \Omega_2$
    max_temp = 100; t_rise = 100
    ch = ConstraintHandler(dh);
    ∂Ω₁ = union(getfaceset.((grid,), ["left", "right"])...)
    add!(ch, Dirichlet(:u, ∂Ω₁, (x, t) -> 0));
    ∂Ω₂ = union(getfaceset.((grid,), ["top", "bottom"])...)
    add!(ch, Dirichlet(:u, ∂Ω₂, (x, t) -> max_temp * clamp(t / t_rise, 0, 1)))
    close!(ch)

    # Create and return the `FEDefinition`
    return FEDefinition(;dh=dh, ch=ch, cv=cellvalues, m=FicksLaw())
end;

struct PostProcessing{PVD}
    pvd::PVD
end
PostProcessing() = PostProcessing(paraview_collection("transient-heat.pvd"));

function FESolvers.postprocess!(post::PostProcessing, p, step, solver)
    @info "postprocessing step $step"
    dh = FP.getdh(p)
    vtk_grid("transient-heat-$step", dh) do vtk
        vtk_point_data(vtk, dh, FP.getunknowns(p))
        vtk_save(vtk)
        post.pvd[step] = vtk
    end
end;

function FESolvers.close_problem(post::PostProcessing, p)
    vtk_save(post.pvd)
end;

def = create_definition()
post = PostProcessing()
problem = FerriteProblem(def, post)
solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper(collect(0.0:1.0:200)));

solve_problem!(solver, problem);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

