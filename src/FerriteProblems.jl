module FerriteProblems

using Printf
using FileIO, JLD2
using Ferrite
using FESolvers, FerriteAssembly, FerriteNeumann
using MaterialModelsBase
import FESolvers: getunknowns, getresidual, getjacobian 
import FerriteAssembly: get_material
export FerriteProblem, FEDefinition, FerriteIO

include("FerritePRs/include_prs.jl")
include("ConvergenceCriteria.jl")
include("FEDefinition.jl")
include("FEBuffer.jl")
include("IO.jl")

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
for op = (:getdh, :getch, :getnh)
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
    :getneumannforce,
    :getmaterial, :getassemblybuffer,
    :get_tolerance_scaling,
    :gettime, :getoldtime, :update_time!, :settime!, 
    :getstate, :getoldstate)
    eval(quote
        $op(p::FerriteProblem, args...) = $op(p.buf, args...)
    end)
end

function FerriteAssembly.update_states!(p::FerriteProblem)
    update_states!(getoldstate(p), getstate(p))
end

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
    
    K = FESolvers.getjacobian(p)
    r = FESolvers.getresidual(p)
    scaling = get_tolerance_scaling(p).assemscaling
    assembler = KeReAssembler(K, r; ch=getch(p), apply_zero=true, scaling=scaling)
    doassemble!(assembler, getstate(p), getassemblybuffer(p); 
        a=a, aold=getoldunknowns(p), old_states=getoldstate(p), Δt=gettime(p)-getoldtime(p)
        )
    r .-= getneumannforce(p)
    apply_zero!(K, r, getch(p))
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