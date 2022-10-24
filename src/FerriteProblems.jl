module FerriteProblems

using Printf
using FileIO, JLD2
using Ferrite
using FESolvers, FerriteAssembly, FerriteNeumann
using MaterialModelsBase
import FESolvers: ScalarWrapper

export FerriteProblem, FEDefinition, safesolve
export savedofdata!, savenodedata!, savecelldata!, saveipdata!, saveglobaldata!
export getdofdata, getnodedata, getcelldata, getipdata, getglobaldata
export initial_conditions!

include("initial_conditions.jl")    # Could maybe add tests and PR into Ferrite - but no time for that now...
include("FEDefinition.jl")
include("FEBuffer.jl")
include("IO.jl")
include("MaterialModelsBase.jl")

struct FerriteProblem{POST,DEF<:FEDefinition,BUF<:FEBuffer,IOT}
    def::DEF
    post::POST
    buf::BUF
    io::IOT
end

function FerriteProblem(def::FEDefinition, post=nothing, io=nothing)
    buf = FEBuffer(def)
    FerriteProblem(def,post,buf,io)
end
function FerriteProblem(def::FEDefinition, post, savefolder::String)
    FerriteProblem(def, post, FerriteIO(savefolder, def, post))
end

"""
    safesolve(solver, def::FEDefinition, args...; kwargs...)

Starts by creating `FerriteProblem` from `def` and `args/kwargs` 
and then wraps the call to `FESolvers.solve_problem!` in 
`try ... finally` to ensure that `close_problem` is called even 
in the case of no convergence.  
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

# FEDefinition: Make functions work directly on `problem`:
for op = (:getdh, :getch, :getnh, :getcv, :getmaterial)
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

# FerriteIO: Make close_problem work on the problem
function close_problem(p::FerriteProblem)
    return close_problem(p.io, p.post)
end
close_problem(::Nothing, args...) = nothing # Do nothing if there is no io defined. 

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