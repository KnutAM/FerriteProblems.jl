struct FEBuffer{T,KT<:AbstractMatrix{T},CC,ST}
    K::KT
    x::Vector{T}
    r::Vector{T}
    f::Vector{T}
    xold::Vector{T}
    cellcache::CC
    state::ST
    old_state::ST
    time::ScalarWrapper{T}
    old_time::ScalarWrapper{T}
end

function FEBuffer(def::FEDefinition)
    n = ndofs(getdh(def))
    K = create_sparsity_pattern(getdh(def))
    x, r, f, xold = ntuple(_->zeros(n), 4)
    cellcache = CellCache(getdh(def))
    state = deepcopy(def.initialstate)
    old_state = deepcopy(def.initialstate)
    time = ScalarWrapper(0.0)
    old_time = ScalarWrapper(0.0)
    return FEBuffer(K, x, r, f, xold, cellcache, state, old_state, time, old_time)
end

# Standard get functions
FESolvers.getjacobian(b::FEBuffer) = b.K
FESolvers.getunknowns(b::FEBuffer) = b.x
FESolvers.getresidual(b::FEBuffer) = b.r
getneumannforce(b::FEBuffer) = b.f 
getcellcache(b::FEBuffer) = b.cellcache

# Variables that will also be updated via special functions
# Unknowns 
getoldunknowns(b::FEBuffer) = b.xold
update_unknowns!(b::FEBuffer) = copy!(b.xold, b.x)

# State variables
getstate(b::FEBuffer) = b.state
getoldstate(b::FEBuffer) = b.old_state

copy_states!(to::Tuple, from::Tuple) = copy_states!.(to, from)
function copy_states!(to::T, from::T) where T<:Vector{<:Vector}
    for (toval, fromval) in zip(to, from)
        copy!(toval, fromval)
    end
end
update_states!(b::FEBuffer) = copy_states!(b.old_state, getstate(b))
reset_states!(b::FEBuffer) = copy_states!(b.state, getoldstate(b))

# Time
gettime(b::FEBuffer) = b.time[]
getoldtime(b::FEBuffer) = b.old_time[]
settime!(b::FEBuffer, new_time) = (b.time[] = new_time)
update_time!(b::FEBuffer) = (b.old_time[] = gettime(b))
