using Ferrite, FerriteProblems, FerriteAssembly, FESolvers
import FerriteProblems as FP
import FerriteAssembly as FA
import FerriteAssembly.ExampleElements: TransientFourier

material() = TransientFourier(#=k=#1.0-3, #=c=#1.0)

function create_definition()
    grid = generate_grid(Quadrilateral, (100, 100));

    cellvalues = CellScalarValues(
        QuadratureRule{2, RefCube}(2),
        Lagrange{2, RefCube, 1}());

    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)

    # Boundary conditions
    # Zero pressure on ∂Ω₁ and linear ramp followed by constant pressure on ∂Ω₂
    max_temp = 100; t_rise = 100
    ch = ConstraintHandler(dh);
    ∂Ω₁ = union(getfaceset.((grid,), ["left", "right"])...)
    add!(ch, Dirichlet(:u, ∂Ω₁, (x, t) -> 0));
    ∂Ω₂ = union(getfaceset.((grid,), ["top", "bottom"])...)
    add!(ch, Dirichlet(:u, ∂Ω₂, (x, t) -> max_temp * clamp(t / t_rise, 0, 1)))
    close!(ch)

    # Body load: constant heat source, f=5.0e-1
    lh = LoadHandler(dh)
    add!(lh, BodyLoad(:u, 2, Returns(5.0e-1)))

    # Create and return the `FEDefinition`
    domainspec = DomainSpec(dh, material(), cellvalues)
    return FEDefinition(domainspec; ch, lh)
end;

struct TH_PostProcessing{PVD}
    pvd::PVD
end
TH_PostProcessing() = TH_PostProcessing(paraview_collection("transient-heat.pvd"));

function FESolvers.postprocess!(post::TH_PostProcessing, p, step, solver)
    if step < 5 || mod(step, 20) == 0
        @info "postprocessing step $step"
        dh = FP.get_dofhandler(p)
        vtk_grid("transient-heat-$step", dh) do vtk
            vtk_point_data(vtk, dh, FP.getunknowns(p))
            vtk_save(vtk)
            post.pvd[step] = vtk
        end
    end
end;

function FP.close_postprocessing(post::TH_PostProcessing, p)
    vtk_save(post.pvd)
end;

def = create_definition()
post = TH_PostProcessing()
problem = FerriteProblem(def, post)
solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper(collect(0.0:1.0:200)));

solve_problem!(problem, solver);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
