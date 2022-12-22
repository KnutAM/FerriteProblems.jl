module FerriteProblems

using Printf
using FileIO, JLD2
using Ferrite
using FESolvers, FerriteAssembly, FerriteNeumann
using MaterialModelsBase
import FESolvers: ScalarWrapper
import FESolvers: getunknowns, getresidual, getjacobian 

export FerriteProblem, FEDefinition, FerriteIO

include("FerritePRs/include_prs.jl")
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
    FESolvers.postprocess!(p::FerriteProblem, step, solver)

When `FESolvers` call this function for `p::FerriteProblem`, 
the following function

    FESolvers.postprocess!(post, p::FerriteProblem, step, solver)

is called where `post=p.post` (unless you define a different override). 
This allows you to easily define the dispatch on your postprocessing 
type as `FESolvers.postprocess!(post::MyPostType, p, step, solver)`
Note that the `solver` input argument is required, but can be 
accounted for by defining, e.g. 
`FESolvers.postprocess!(post::MyPostType, p, step, args...)`
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

"""
    close_postprocessing(post::MyPostType, p::FerriteProblem)

This function is called to close any open files manually created during 
the postprocessing with the custom postprocessing type `MyPostType`. 
Note that the file streams in p.io::FerriteIO are 
automatically closed and don't require any special handling.
"""
close_postprocessing(args...) = nothing

# FerriteIO: Make close_problem work on the problem, see IO.jl for definition
function FESolvers.close_problem(p::FerriteProblem)
    close_postprocessing(p.post, p)  # User defined closing for custom postprocessing
    close_io(p.io, p.post)           # FerriteIO function for closing those streams
end

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