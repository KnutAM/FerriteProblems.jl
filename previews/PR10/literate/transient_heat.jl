# # Linear Time Dependent Problem
# This example is the same example as 
# [`FESolvers.jl`'s transient heat flow](https://knutam.github.io/FESolvers.jl/dev/examples/transient_heat/) 
# which was taken from 
# [`Ferrite.jl`'s transient heat flow](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/transient_heat_equation/).
# Please see the theoretical derivations in those examples, with the specific formulation used here in the former. 
# 
# ## Commented Program
# 
# Now we solve the problem by using FerriteProblems. 
#md # The full program, without comments, can be found in the next [section](@ref transient_heat_equation-plain-program).
# 
# First we load required packages
using Ferrite, FerriteProblems, FerriteAssembly, FESolvers
import FerriteProblems as FP
import FerriteAssembly as FA
import FerriteAssembly.ExampleElements: TransientFourier

# ## Physics
# The transient heat element, `TransientFourier`, is available 
# from `FerriteAssembly.ExampleElements`, noting that in the 
# original example, the heat capacity, c=1, is implicitly assumed.
material() = TransientFourier(#=k=#1.0-3, #=c=#1.0)

# The material doesn't include the heat source: 
# We later add this with `FerriteAssembly`'s 
# `LoadHandler`, which works similar 
# to the ConstraintHandler.

# ## Problem setup
# We start by a function that will create the problem definition
function create_definition()
    grid = generate_grid(Quadrilateral, (100, 100));

    cellvalues = CellScalarValues(
        QuadratureRule{2, RefCube}(2), 
        Lagrange{2, RefCube, 1}());

    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
    
    ## Boundary conditions
    ## Zero pressure on ∂Ω₁ and linear ramp followed by constant pressure on ∂Ω₂
    max_temp = 100; t_rise = 100
    ch = ConstraintHandler(dh);
    ∂Ω₁ = union(getfaceset.((grid,), ["left", "right"])...)
    add!(ch, Dirichlet(:u, ∂Ω₁, (x, t) -> 0));
    ∂Ω₂ = union(getfaceset.((grid,), ["top", "bottom"])...)
    add!(ch, Dirichlet(:u, ∂Ω₂, (x, t) -> max_temp * clamp(t / t_rise, 0, 1)))
    close!(ch)

    ## Body load: constant heat source, f=5.0e-1
    lh = LoadHandler(dh)
    add!(lh, BodyLoad(:u, 2, Returns(5.0e-1)))
    
    ## Create and return the `FEDefinition`
    domainspec = DomainSpec(dh, material(), cellvalues)
    return FEDefinition(domainspec; ch, lh)
end;

# ## Postprocessing
# After defining all the physics and problem setup, we must decide what data to save.
# In this example, we use the vtk-file exports as in the original example. To this end, 
# we define the custom postprocessing struct
struct TH_PostProcessing{PVD}
    pvd::PVD
end
TH_PostProcessing() = TH_PostProcessing(paraview_collection("transient-heat.pvd"));

# And the postprocessing function that is called after each time step
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

# At the end of the simulation, we want to finish all IO operations. 
# We can then define the function `close_postprocessing` which will be called 
# even in the case that an error is thrown during the simulation
function FP.close_postprocessing(post::TH_PostProcessing, p)
    vtk_save(post.pvd)
end;

# And now we create the problem type, and define the QuasiStaticSolver with 
# the LinearProblemSolver as well as fixed time steps 
def = create_definition()
post = TH_PostProcessing()
problem = FerriteProblem(def, post)
solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper(collect(0.0:1.0:200)));

# Finally, we can solve the problem
solve_problem!(problem, solver);

#md # ## [Plain program](@id transient_heat_equation-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here:
#md # [`transient_heat.jl`](transient_heat.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
