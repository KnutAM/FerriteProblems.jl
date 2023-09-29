module FerriteProblems
import Base: @kwdef # Support julia <1.9import Base: @kwdef # Support julia <1.9
using Printf
using FileIO, JLD2
using Ferrite
using FESolvers, FerriteAssembly
import FESolvers: getunknowns, getresidual, getjacobian 
export FerriteProblem, FEDefinition, FerriteIO
import FerriteAssembly: get_dofhandler, get_material, get_state, get_old_state

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
    FerriteProblem(def::FEDefinition, post=nothing, io=nothing; kwargs...)
    FerriteProblem(def::FEDefinition, post, savefolder::String; kwargs...)
    
Create a FerriteProblem from [`def`](@ref FEDefinition).
Postprocessing can be added as `post`, see [`FESolvers.postprocess!`](@ref).
File input/output using [`FerriteIO`](@ref) can be added with `io`. 
It is also possible to just give the `savefolder`, i.e. where 
to save the output when the default `FerriteIO` can be used.

Supported keyword arguments are
* `autodiffbuffer::Bool`: Should `FerriteAssembly.jl`'s `AutoDiffCellBuffer` be used? 
  This will make the assembly faster if automatic differentiation is used, and can also be used 
  without automatic differentiation (but with a slight extra computational overhead)
* `threading::Bool`: Should threading be used? 
"""
function FerriteProblem(def::FEDefinition, post=nothing, io=nothing; kwargs...)
    buf = FEBuffer(def; kwargs...)
    FerriteProblem(def, post, buf, io)
end
function FerriteProblem(def::FEDefinition, post, savefolder::String; kwargs...)
    FerriteProblem(def, post, FerriteIO(savefolder, def, post); kwargs...)
end

function Base.show(io::IO, p::FerriteProblem)
    n_dofs = ndofs(get_dofhandler(p))
    print(io, "FerriteProblem[ndofs=$n_dofs]")
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
function FESolvers.postprocess!(::Nothing, ::FerriteProblem, args...)
    return nothing
end

# FEDefinition: Make functions work directly on `problem`:
for op = (:get_dofhandler, :get_constrainthandler, :get_loadhandler)
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
    :get_old_unknowns, :update_unknowns!, 
    :get_external_force,
    :get_material, :getassemblybuffer,
    :get_tolerance_scaling,
    :get_time, :get_old_time, :update_time!, :set_time!, 
    :get_state, :get_old_state)
    eval(quote
        $op(p::FerriteProblem, args...) = $op(p.buf, args...)
    end)
end

function FerriteAssembly.update_states!(p::FerriteProblem)
    update_states!(getassemblybuffer(p))
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

addstep!(io::FerriteIO, p::FerriteProblem) = addstep!(io, get_time(p))

# Actual FE-"work"
# * update_to_next_step!
# * update_problem!
# * calculate_convergence_measure
# * handle_converged!

function FESolvers.update_to_next_step!(p::FerriteProblem, time)
    ch = get_constrainthandler(p)
    # Update the current time
    set_time!(p, time)
    set_time_increment!(getassemblybuffer(p), get_time(p)-get_old_time(p))

    # Apply constraints, including Dirichlet BC
    update!(ch, time)
    apply!(FESolvers.getunknowns(p), ch)
    
    # Apply external load
    f = get_external_force(p)
    fill!(f, 0)
    apply!(f, get_loadhandler(p), time)
    apply_zero!(f, ch) # Make force zero at constrained dofs (to be compatible with apply local)
end

function FESolvers.update_problem!(p::FerriteProblem, Δa, update_spec)
    # Update a if Δa is given
    a = FESolvers.getunknowns(p)
    ch = get_constrainthandler(p)
    if !isnothing(Δa)
        apply_zero!(Δa, ch)
        a .+= Δa
    end
    
    # If `update_spec` requests, update the `m=set_jacobian_type(m, type)`
    # This defaults to no-op if not defined. 
    change_material_jacobian_type!(p.buf, FESolvers.get_update_type(update_spec))

    K = FESolvers.getjacobian(p)
    r = FESolvers.getresidual(p)
    if FESolvers.should_update_jacobian(update_spec) # Update both residual and jacobian
        scaling = get_tolerance_scaling(p).assemscaling
        assembler = KeReAssembler(K, r; scaling=scaling)
    elseif FESolvers.should_update_residual(update_spec) # update only residual
        assembler = ReAssembler(r; scaling=get_tolerance_scaling(p).assemscaling)
    else # Only update a
        return nothing
    end
    map!(-, r, get_external_force(p))
    work!(assembler, getassemblybuffer(p); a, aold=get_old_unknowns(p))
    apply_zero!(K, r, ch)
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