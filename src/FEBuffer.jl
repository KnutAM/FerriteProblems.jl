"""
    FEBuffer(K,x,r,f,xold,cellbuffer,state,old_state,time,old_time)

A buffer to hold all values that are required to simulate, 
but that are uniqely defined from the simulation definition
"""
mutable struct FEBuffer{T,KT<:AbstractMatrix{T},CB,ST,TS} # internal 
    const K::KT
    const x::Vector{T}
    const r::Vector{T}
    const f::Vector{T}
    const xold::Vector{T}
    const cellbuffer::CB
    const state::ST
    const old_state::ST
    time::T
    old_time::T
    const tolscaling::TS
end

"""
    cellbuffertype(def)

By default, `cellbuffertype(::Any) = CellBuffer`. 
However, if you are using the automatic differentiation, 
much better assembly speed can be achieved by defining
```
FerriteProblems.cellbuffertype(::FEDefinition) = AutoDiffCellBuffer
```
In this case, any defined `element_routine!` should 
use `getCellBuffer(buffer)` to get a `CellBuffer` instead 
of `AutoDiffCellBuffer`. 
This is not necessary in `element_residual!`
"""
cellbuffertype(::Any) = CellBuffer

# makecellbuffer could also be used to use a custom `AbstractCellBuffer`
makecellbuffer(def) = makecellbuffer(def, dothreaded(def), cellbuffertype(def))
makecellbuffer(def, threaded::Val{false}, CB::Type{<:CellBuffer}) = setup_cellbuffer(getdh(def), getcv(def), getmaterial(def), getbodyload(def), allocate_material_cache(def))
makecellbuffer(def, threaded::Val{false}, CB::Type{<:AutoDiffCellBuffer}) = setup_ad_cellbuffer(def.initialstate, getdh(def), getcv(def), getmaterial(def), getbodyload(def), allocate_material_cache(def))
makecellbuffer(def, threaded::Val{true}, CB) = create_threaded_CellBuffers(makecellbuffer(def, Val{false}(), CB))

function FEBuffer(def::FEDefinition)
    n = ndofs(getdh(def))
    K = create_sparsity_pattern(getdh(def))
    x, r, f = (zeros(n) for _ in 1:3)
    foreach(ic->apply_analytical!(x, getdh(def), ic[1], ic[2]), pairs(def.ic))
    xold = deepcopy(x)
    cellbuffer = makecellbuffer(def)
    state = deepcopy(def.initialstate)
    old_state = deepcopy(def.initialstate)
    time = 0.0
    old_time = 0.0
    return FEBuffer(K, x, r, f, xold, cellbuffer, state, old_state, time, old_time, TolScaling(def.cc, def))
end

# Standard get functions
"""
    FerriteProblems.getjacobian(p::FerriteProblem)

Get the current jacobian matrix from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FerriteProblems.getjacobian`
"""
getjacobian(b::FEBuffer) = b.K

"""
    FerriteProblems.getunknowns(p::FerriteProblem)

Get the current vector of unknowns from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FerriteProblems.getunknowns`
"""
getunknowns(b::FEBuffer) = b.x

"""
    FerriteProblems.getresidual(p::FerriteProblem)

Get the current residual vector from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FerriteProblems.getresidual`
"""
getresidual(b::FEBuffer) = b.r

"""
    FerriteProblems.getneumannforce(p::FerriteProblem)

Get the current external force vector caused by 
Neumann boundary conditions. Note that this vector 
does not include external forces added during the 
cell assembly; only forces added with the `NeumannHandler`
"""
getneumannforce(b::FEBuffer) = b.f 

"""
    FerriteProblems.getcellbuffer(p::FerriteProblem)

Get the cell buffers used during the assembly.
"""
getcellbuffer(b::FEBuffer) = b.cellbuffer   # Internal

"""
    FerriteProblems.get_tolerance_scaling(p::FerriteProblem)

Get the `TolScaling` type that controls the convergence 
measure to be compared with the solver's tolerance
"""
get_tolerance_scaling(b::FEBuffer) = b.tolscaling


# Variables that will also be updated via special functions
# Unknowns 
"""
    FerriteProblems.getoldunknowns(p::FerriteProblem)

Get the vector of unknowns from the previously converged step
"""
getoldunknowns(b::FEBuffer) = b.xold 

"""
    FerriteProblems.update_unknowns!(p::FerriteProblem)

Update the vector of "old" unknowns to the values of the current vector of unknowns
"""
update_unknowns!(b::FEBuffer) = copy!(b.xold, b.x) # internal

# State variables
"""
    FerriteProblems.getstate(p::FerriteProblem)

Get the current state variables
"""
getstate(b::FEBuffer) = b.state

"""
    FerriteProblems.getoldstate(p::FerriteProblem)

Get the state variables from the previously converged step
"""
getoldstate(b::FEBuffer) = b.old_state

"""
    FerriteProblems.copy_states!(to, from)

Perform a "`deepcopy!`" of states in `from` into `to`.
"""
copy_states!(to::Tuple, from::Tuple) = copy_states!.(to, from)  # internal
function copy_states!(to::T, from::T) where T<:Vector{<:Vector}
    for (toval, fromval) in zip(to, from)
        copy_states!(toval, fromval)
    end
end
function copy_states!(to::T, from::T) where T<:Dict{<:Union{Int,String},<:Vector}
    for (key, fromval) in from
        copy_states!(to[key], fromval)
    end
end
copy_states!(::T, ::T) where T<:Union{Vector{Nothing},Dict{Int,Nothing}} = nothing 

@inline copy_states!(to::T, from::T) where T<:Vector{ET} where ET = copy_states!(Val{isbitstype(ET)}(), to, from)
@inline copy_states!(::Val{true}, to::Vector, from::Vector) = copy!(to,from)
@inline copy_states!(::Val{false}, to::Vector, from::Vector) = copy!(to,deepcopy(from))

"""
    FerriteProblems.update_states!(p::FerriteProblem)

Update the "old" state variables to the current values.
Called after convergence
"""
update_states!(b::FEBuffer) = copy_states!(b.old_state, getstate(b)) # internal

"""
    FerriteProblems.reset_states!(p::FerriteProblem)

Reset the current state variables to the old state values. 
Called after each solution iteration unless the solution has converged. 
"""
reset_states!(b::FEBuffer) = copy_states!(b.state, getoldstate(b))  # internal

# Time
"""
    FerriteProblems.gettime(p::FerriteProblem)

Get the current time
"""
gettime(b::FEBuffer) = b.time

"""
    FerriteProblems.getoldtime(p::FerriteProblem)

Get time of the previous converged step
"""
getoldtime(b::FEBuffer) = b.old_time

"""
    FerriteProblems.settime!(p::FerriteProblem, new_time)

Set the current time to `new_time`
Called when starting a new step (or when attempting the same 
step number with a new time increment)
"""
settime!(b::FEBuffer, new_time) = (b.time = new_time) # internal

"""
    FerriteProblems.update_time!(p::FerriteProblem)

Update the old time to the current time. 
Called after convergence
"""
update_time!(b::FEBuffer) = (b.old_time = gettime(b)) # internal