module FerriteProblems

using Printf
using FileIO, JLD2
using Ferrite
using FESolvers, FerriteAssembly, FerriteNeumann
using MaterialModelsBase
import FESolvers: ScalarWrapper
import FESolvers: getunknowns, getresidual, getjacobian 

const FP = FerriteProblems # Provide shorthand instead of exporting all functions
export FP
export FerriteProblem, FEDefinition, FerriteIO
export safesolve
export initial_conditions!

include("FerritePRs.jl")
include("initial_conditions.jl")    # Should be updated and added to FerritePRs.jl...
include("ConvergenceCriteria.jl")
include("FEDefinition.jl")
include("FEBuffer.jl")
include("IO.jl")
include("MaterialModelsBase.jl")

"""
    FerriteProblem(def::FEDefinition, post, buf::FEBuffer, io::[FerriteIO])
    
The main problem type that holds all variables to solve a particular problem 
using `FESolvers`
"""
struct FerriteProblem{POST,DEF<:FEDefinition,BUF<:FEBuffer,IOT}
    def::DEF
    post::POST
    buf::BUF
    io::IOT
end

"""
    FerriteProblem(def::FEDefinition, post=nothing, io=nothing)
    
Constructor that makes the minimum required to run simulation. 
Optional postprocessing and `io::FerriteIO` if desired. 
"""
function FerriteProblem(def::FEDefinition, post=nothing, io=nothing)
    buf = FEBuffer(def)
    FerriteProblem(def,post,buf,io)
end

"""
    FerriteProblem(def::FEDefinition, post, savefolder::String)

Constructor with automatic generation of `FerriteIO` by just specifying 
in which folder data should be saved. 
"""
function FerriteProblem(def::FEDefinition, post, savefolder::String)
    FerriteProblem(def, post, FerriteIO(savefolder, def, post))
end

function Base.show(io::IO, p::FerriteProblem)
    print(io, "FerriteProblem[ndofs=$(ndofs(getdh(p)))]")
end

"""
    safesolve(solver, def::FEDefinition, args...; kwargs...)

Starts by creating `FerriteProblem` from `def` and `args/kwargs` 
and then wraps the call to `FESolvers.solve_problem!` in 
`try ... finally` to ensure that [`close_problem`](@ref) is 
called even in the case of no convergence.  
"""
function safesolve(solver, def::FEDefinition, args...; kwargs...)
    problem = FerriteProblem(def, args...; kwargs...)
    try
        solve_problem!(solver, problem)
    finally
        close_problem(problem)
    end
    return problem
end

"""
    FESolvers.postprocess!(p::FerriteProblem, step, solver)

When `FESolvers` call this function for `p::FerriteProblem`, 
the following function

    FESolvers.postprocess!(post, p::FerriteProblem, step, solver)

is called where `post=p.post` (unless you define a different override). 
This allows you to easily define the dispatch on your postprocessing 
type as `FESolvers.postprocess!(post::MyPostType, p, step, solver)`
"""
function FESolvers.postprocess!(p::FerriteProblem, step, solver)
    return FESolvers.postprocess!(p.post, p, step, solver)
end

# FEDefinition: Make functions work directly on `problem`:
for op = (:getdh, :getch, :getnh, :getcv, :getmaterial, :getbodyload, :dothreaded)
    eval(quote
        $op(p::FerriteProblem) = $op(p.def)
    end)
end

# FEBuffer: Make functions work directly on `problem`:
for op = (:getunknowns, :getresidual, :getjacobian)
    eval(quote
        FESolvers.$op(p::FerriteProblem) = FESolvers.$op(p.buf)
    end)
end

for op = (
    :getoldunknowns, :update_unknowns!, 
    :getneumannforce, :getcellbuffer, 
    :get_tolerance_scaling,
    :gettime, :getoldtime, :update_time!,
    :getstate, :getoldstate, :update_states!, :reset_states!)
    eval(quote
        $op(p::FerriteProblem) = $op(p.buf)
    end)
end
settime!(p::FerriteProblem, t) = settime!(p.buf, t)

# FerriteIO: Make close_problem work on the problem, see IO.jl for definition
"""
    close_problem(p::FerriteProblem)

Method for closing all open files before ending the simulation
"""
function FESolvers.close_problem(p::FerriteProblem)
    FESolvers.close_problem(p.post, p)  # User defined closing for custom postprocessing
    return close_problem(p.io, p.post)
end
close_problem(::Any, args...) = nothing # Do nothing if there is no io defined. 

addstep!(io::FerriteIO, p::FerriteProblem) = addstep!(io, gettime(p))

# Actual FE-"work"
# * update_to_next_step!
# * update_problem!
# * calculate_convergence_measure
# * handle_converged!

function FESolvers.update_to_next_step!(p::FerriteProblem, time)
    # Update the current time
    settime!(p, time)
    
    # Apply Neumann BC
    f = getneumannforce(p)
    fill!(f, 0)
    apply!(f, getnh(p), time)
    
    # Apply constraints, including Dirichlet BC
    update!(getch(p), time)
    apply!(FESolvers.getunknowns(p), getch(p))
end

function FESolvers.update_problem!(p::FerriteProblem, Δa; kwargs...)
    # Update a if Δa is given
    a = FESolvers.getunknowns(p)
    if !isnothing(Δa)
        apply_zero!(Δa, getch(p))
        a .+= Δa
    end
    reset_states!(p)    # "state = old_state"
    state = getstate(p)
    Δt = gettime(p) - getoldtime(p)
    K = FESolvers.getjacobian(p)
    r = FESolvers.getresidual(p)
    aold = getoldunknowns(p)
    scaling = get_tolerance_scaling(p).assemscaling
    _doassemble!(K, r, getcellbuffer(p), state, getdh(p), a, aold, Δt, scaling, p, dothreaded(p))
    r .-= getneumannforce(p)
    apply_zero!(K, r, getch(p))
end

function _doassemble!(K, r, cellbuffer, state, dh, a, aold, Δt, scaling, p, threaded::Val{false})
    assembler = start_assemble(K, r)
    reset_scaling!(scaling)
    FerriteAssembly.doassemble!(assembler, cellbuffer, state, dh, a, aold, Δt, scaling)
end

function _doassemble!(K, r, cellbuffers, states, dh, a, aold, Δt, scalings, p, threaded::Val{true})
    colors = p.def.colors
    assemblers = create_threaded_assemblers(K, r);
    reset_scaling!.(scalings)
    doassemble!(assemblers, cellbuffers, states, dh, colors, a, aold, Δt, scalings)
end


function FESolvers.calculate_convergence_measure(p::FerriteProblem, Δa, iter)
    ts = get_tolerance_scaling(p)
    r = FESolvers.getresidual(p)
    return FESolvers.calculate_convergence_measure(ts, r, Δa, iter, p)
end

function FESolvers.handle_converged!(p::FerriteProblem)
    # Set old = current
    update_time!(p)
    update_states!(p)
    update_unknowns!(p)
end

end