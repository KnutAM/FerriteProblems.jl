"""
    FEBuffer(K,x,r,f,xold,cellbuffer,state,old_state,time,old_time)

A buffer to hold all values that are required to simulate, 
but that are uniqely defined from the simulation definition
"""
struct FEBuffer{T,KT<:AbstractMatrix{T},CB,ST} # internal 
    K::KT
    x::Vector{T}
    r::Vector{T}
    f::Vector{T}
    xold::Vector{T}
    cellbuffer::CB
    state::ST
    old_state::ST
    time::ScalarWrapper{T}
    old_time::ScalarWrapper{T}
end

function FEBuffer(def::FEDefinition)
    n = ndofs(getdh(def))
    K = create_sparsity_pattern(getdh(def))
    x, r, f = (zeros(n) for _ in 1:3)
    foreach(ic->initial_conditions!(x, getdh(def), ic[1], ic[2]), pairs(def.ic))
    xold = deepcopy(x)
    cellbuffer = CellBuffer(getdh(def), getcv(def), getmaterial(def), getbodyload(def), allocate_material_cache(def))
    state = deepcopy(def.initialstate)
    old_state = deepcopy(def.initialstate)
    time = ScalarWrapper(0.0)
    old_time = ScalarWrapper(0.0)
    return FEBuffer(K, x, r, f, xold, cellbuffer, state, old_state, time, old_time)
end

# Standard get functions
"""
    FP.getjacobian(p::FerriteProblem)

Get the current jacobian matrix from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FP.getjacobian`
"""
getjacobian(b::FEBuffer) = b.K

"""
    FP.getunknowns(p::FerriteProblem)

Get the current vector of unknowns from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FP.getunknowns`
"""
getunknowns(b::FEBuffer) = b.x

"""
    FP.getresidual(p::FerriteProblem)

Get the current residual vector from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FP.getresidual`
"""
getresidual(b::FEBuffer) = b.r

"""
    FP.getneumannforce(p::FerriteProblem)

Get the current external force vector caused by 
Neumann boundary conditions. Note that this vector 
does not include external forces added during the 
cell assembly; only forces added with the `NeumannHandler`
"""
getneumannforce(b::FEBuffer) = b.f 

"""
    FP.getcellbuffer(p::FerriteProblem)

Get the cell buffers used during the assembly.
"""
getcellbuffer(b::FEBuffer) = b.cellbuffer   # Internal

# Variables that will also be updated via special functions
# Unknowns 
"""
    FP.getoldunknowns(p::FerriteProblem)

Get the vector of unknowns from the previously converged step
"""
getoldunknowns(b::FEBuffer) = b.xold 

"""
    FP.update_unknowns!(p::FerriteProblem)

Update the vector of "old" unknowns to the values of the current vector of unknowns
"""
update_unknowns!(b::FEBuffer) = copy!(b.xold, b.x) # internal

# State variables
"""
    FP.getstate(p::FerriteProblem)

Get the current state variables
"""
getstate(b::FEBuffer) = b.state

"""
    FP.getoldstate(p::FerriteProblem)

Get the state variables from the previously converged step
"""
getoldstate(b::FEBuffer) = b.old_state

"""
    FP.copy_states!(to, from)

Perform a "`deepcopy!`" of states in `from` into `to`.
"""
copy_states!(to::Tuple, from::Tuple) = copy_states!.(to, from)  # internal
function copy_states!(to::T, from::T) where T<:Vector{<:Vector}
    for (toval, fromval) in zip(to, from)
        copy_states!(toval, fromval)
    end
end
function copy_states!(to::T, from::T) where T<:Dict{Int,<:Vector}
    for (key, fromval) in from
        copy_states!(to[key], fromval)
    end
end
@inline copy_states!(to::T, from::T) where T<:Vector{ET} where ET = copy_states!(Val{isbitstype(ET)}(), to, from)
@inline copy_states!(::Val{true}, to::Vector, from::Vector) = copy!(to,from)
@inline copy_states!(::Val{false}, to::Vector, from::Vector) = copy!(to,deepcopy(from))

"""
    FP.update_states!(p::FerriteProblem)

Update the "old" state variables to the current values.
Called after convergence
"""
update_states!(b::FEBuffer) = copy_states!(b.old_state, getstate(b)) # internal

"""
    FP.reset_states!(p::FerriteProblem)

Reset the current state variables to the old state values. 
Called after each solution iteration unless the solution has converged. 
"""
reset_states!(b::FEBuffer) = copy_states!(b.state, getoldstate(b))  # internal

# Time
"""
    FP.gettime(p::FerriteProblem)

Get the current time
"""
gettime(b::FEBuffer) = b.time[]

"""
    FP.getoldtime(p::FerriteProblem)

Get time of the previous converged step
"""
getoldtime(b::FEBuffer) = b.old_time[]

"""
    FP.settime!(p::FerriteProblem, new_time)

Set the current time to `new_time`
Called when starting a new step (or when attempting the same 
step number with a new time increment)
"""
settime!(b::FEBuffer, new_time) = (b.time[] = new_time) # internal

"""
    FP.update_time!(p::FerriteProblem)

Update the old time to the current time. 
Called after convergence
"""
update_time!(b::FEBuffer) = (b.old_time[] = gettime(b)) # internal