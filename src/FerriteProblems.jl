module FerriteProblems
import Base: @kwdef # Support julia <1.9
using Printf
using FileIO, JLD2
using Ferrite
using FESolvers, FerriteAssembly, FerriteNeumann
import FESolvers: getunknowns, getresidual, getjacobian 
import FerriteAssembly: get_material
export FerriteProblem, FEDefinition, FerriteIO

include("FerritePRs/include_prs.jl")
include("ConvergenceCriteria.jl")
include("FEDefinition.jl")
include("FEBuffer.jl")
include("IO.jl")
include("CustomStiffness.jl")

struct FerriteProblem{DEF<:FEDefinition,POST,BUF<:FEBuffer,IOT}
    def::DEF
    post::POST
    buf::BUF
    io::IOT
end

"""
    FerriteProblem(def::FEDefinition, post=nothing, io=nothing)
    
Create a FerriteProblem from [`def`](@ref FEDefinition).
Postprocessing can be added as `post`, see [`FESolvers.postprocess!`](@ref).
File input/output using `FerriteIO` can be added with `io`. 
It is possible to give the folder where to save the output (i.e. io::String), 
or to construct [`FerriteIO`](@ref) with more options directly.
"""
function FerriteProblem(def::FEDefinition, post=nothing, io=nothing)
    buf = FEBuffer(def)
    FerriteProblem(def, post, buf, io)
end
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
    :get_material, :getassemblybuffer,
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
    p.def.fesolverfuns.update_to_next_step!(p, time)
end

function fp_update_to_next_step!(p::FerriteProblem, time)
    # Update the current time
    settime!(p, time)

    # Apply constraints, including Dirichlet BC
    update!(getch(p), time)
    apply!(FESolvers.getunknowns(p), getch(p))
    
    # Apply Neumann BC
    f = getneumannforce(p)
    fill!(f, 0)
    apply!(f, getnh(p), time)
    apply_zero!(f, getch(p)) # Make force zero at constrained dofs (to be compatible with apply local)
end

function FESolvers.update_problem!(p::FerriteProblem, args...; kwargs...)
    p.def.fesolverfuns.update_problem!(p, args...; kwargs...)
end
function fp_update_problem!(p::FerriteProblem, Δa, update_spec)
    # Update a if Δa is given
    a = FESolvers.getunknowns(p)
    if !isnothing(Δa)
        apply_zero!(Δa, getch(p))
        a .+= Δa
    end
    
    # If `update_spec` requests, update the `m=set_jacobian_type(m, type)`
    # This defaults to no-op if not defined. 
    change_material_jacobian_type!(p.buf, FESolvers.get_update_type(update_spec))

    K = FESolvers.getjacobian(p)
    r = FESolvers.getresidual(p)
    if FESolvers.should_update_jacobian(update_spec) # Update both residual and jacobian
        scaling = get_tolerance_scaling(p).assemscaling
        assembler = KeReAssembler(K, r; ch=getch(p), apply_zero=true, scaling=scaling)
    elseif FESolvers.should_update_residual(update_spec) # update only residual
        assembler = ReAssembler(r; scaling=get_tolerance_scaling(p).assemscaling)
    else # Only update a
        return nothing
    end
    map!(-, r, getneumannforce(p))
    doassemble!(assembler, getstate(p), getassemblybuffer(p); 
        a=a, aold=getoldunknowns(p), old_states=getoldstate(p), Δt=gettime(p)-getoldtime(p)
        )
    return nothing
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