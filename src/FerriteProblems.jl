module FerriteProblems

using Printf
using FileIO, JLD2
using Ferrite
using FESolvers, FerriteAssembly, FerriteNeumann
import FESolvers: ScalarWrapper

export FerriteProblem, FEDefinition
export savedofdata!, savenodedata!, savecelldata!, saveipdata!, saveglobaldata!
export getdofdata, getnodedata, getcelldata, getipdata, getglobaldata

include("FEDefinition.jl")
include("FEBuffer.jl")
include("IO.jl")

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
    :getneumannforce, :getcellcache, 
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
    FerriteAssembly.doassemble!(K, r, a, aold, state, getdh(p), getcv(p), getmaterial(p), Δt, getcellcache(p))
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