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

include("initial_conditions.jl")    # Could maybe add tests and PR into Ferrite - but no time for that now...
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
for op = (:getdh, :getch, :getnh, :getcv, :getmaterial, :getbodyload)
    eval(quote
        FerriteProblems.$op(p::FerriteProblem) = FerriteProblems.$op(p.def)
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
function close_problem(p::FerriteProblem)
    return close_problem(p.io, p.post)
end
close_problem(::Nothing, args...) = nothing # Do nothing if there is no io defined. 

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

function FESolvers.update_problem!(p::FerriteProblem, Δa=nothing)
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
    assembler = start_assemble(K, r)
    FerriteAssembly.doassemble!(assembler, getcellbuffer(p), state, getdh(p), a, aold, Δt)
    r .-= getneumannforce(p)
    apply_zero!(K, r, getch(p))
end

function FESolvers.calculate_convergence_measure(p::FerriteProblem)
    r = FESolvers.getresidual(p)
    return sqrt(sum(i->r[i]^2, Ferrite.free_dofs(getch(p))))
end

function FESolvers.handle_converged!(p::FerriteProblem)
    # Set old = current
    update_time!(p)
    update_states!(p)
    update_unknowns!(p)
end

end