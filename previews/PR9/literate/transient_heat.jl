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

# ## Physics
# First, we need to define the material behavior. 
@kwdef struct FourierMaterial{T}
    k::T=1.0e-3    # Thermal conductivity
    f::T=5.0e-1    # Constant heat source
end
# where we could have defined the heat source using the bodyload type 
# available via the cellbuffer, but it is not necessary for a constant 
# heat source. 

# We then define element routine following `FerriteAssembly`
function FerriteAssembly.element_routine!(
    Ke, re, state, ue, m::FourierMaterial, cellvalues, buffer
    )
    Δt = FA.get_time_increment(buffer)
    ue_old = FA.get_aeold(buffer)
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

# ## Problem setup
# We start by a function that will create the problem definition
function create_definition()
    ## **Grid**
    grid = generate_grid(Quadrilateral, (100, 100));

    ## **Cell values**
    cellvalues = CellScalarValues(
        QuadratureRule{2, RefCube}(2), 
        Lagrange{2, RefCube, 1}());

    ## **Degrees of freedom**
    ## After this, we can define the `DofHandler` and distribute the DOFs of the problem.
    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)

    ## **Boundary conditions**
    ## Zero pressure on $\partial \Omega_1$ and linear ramp followed by constant pressure on $\partial \Omega_2$
    max_temp = 100; t_rise = 100
    ch = ConstraintHandler(dh);
    ∂Ω₁ = union(getfaceset.((grid,), ["left", "right"])...)
    add!(ch, Dirichlet(:u, ∂Ω₁, (x, t) -> 0));
    ∂Ω₂ = union(getfaceset.((grid,), ["top", "bottom"])...)
    add!(ch, Dirichlet(:u, ∂Ω₂, (x, t) -> max_temp * clamp(t / t_rise, 0, 1)))
    close!(ch)
    
    ## Create and return the `FEDefinition`
    return FEDefinition(;dh=dh, ch=ch, cellvalues=cellvalues, material=FourierMaterial())
end;

# ## Postprocessing
# After defining all the physics and problem setup, we must decide what data to save.
# In this example, we use the vtk-file exports as in the original example. To this end, 
# we define the custom postprocessing struct
struct PostProcessing{PVD}
    pvd::PVD
end
PostProcessing() = PostProcessing(paraview_collection("transient-heat.pvd"));

# And the postprocessing function that is called after each time step
function FESolvers.postprocess!(post::PostProcessing, p, step, solver)
    if step < 5 || mod(step, 20) == 0
        @info "postprocessing step $step"
        dh = FP.getdh(p)
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
function FP.close_postprocessing(post::PostProcessing, p)
    vtk_save(post.pvd)
end;

# And now we create the problem type, and define the QuasiStaticSolver with 
# the LinearProblemSolver as well as fixed time steps 
def = create_definition()
post = PostProcessing()
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
