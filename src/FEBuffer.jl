"""
    FEBuffer(K,x,r,f,xold,cellbuffer,state,old_state,time,old_time)

A buffer to hold all values that are required to simulate, 
but that are uniqely defined from the simulation definition
"""
mutable struct FEBuffer{T,KT<:AbstractMatrix,VT<:AbstractVector,DB,TS} # internal 
    const K::KT
    const x::VT
    const r::VT
    const f::VT
    const xold::VT
    const domain_buffer::DB # DomainBuffer from FerriteAssembly.setup_domainbuffer[s]
    time::T
    old_time::T
    const tolscaling::TS
end

function FEBuffer(def::FEDefinition; kwargs...)
    dh = get_dofhandler(def)
    ch = get_constrainthandler(def)
    n = ndofs(dh)
    K = create_sparsity_pattern(dh, ch)
    x, r, f = zeros.((n,n,n))
    ic = def.initial_conditions # ::NamedTuple 
    foreach((field_name, fun)->apply_analytical!(x, dh, field_name, fun), keys(ic), values(ic))
    xold = deepcopy(x)
    buffer = _setup_assembly(def.domains; a=x, kwargs...)
    time = 0.0
    old_time = 0.0
    tol_scaling = TolScaling(def.convergence_criterion, def)
    return FEBuffer(K, x, r, f, xold, buffer, time, old_time, tol_scaling)
end

function _setup_assembly(db::DomainSpec; a, kwargs...)
    return setup_domainbuffer(db; a, kwargs...)
end
function _setup_assembly(db::Dict{String,<:DomainSpec}; a, kwargs...)
    return setup_domainbuffers(db; a, kwargs...)
end

# Standard get functions
getassemblybuffer(b::FEBuffer) = b.domain_buffer
"""
    FerriteProblems.get_material(p::FerriteProblem)
    FerriteProblems.get_material(p::FerriteProblem, domain_name::String)

Get the material in `p`. For multiple domains, it is necessary to give the `domain_name`
for where to get the material. Note that this is type-unstable and should be avoided in 
performance-critical code sections. This function belongs to `FerriteAssembly.jl`, 
but can be accessed via `FerriteProblems.get_material`.
"""
get_material(b::FEBuffer, args...) = get_material(getassemblybuffer(b), args...)

"""
    FerriteProblems.getjacobian(p::FerriteProblem)

Get the current jacobian matrix from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FerriteProblems.getjacobian`
"""
FESolvers.getjacobian(b::FEBuffer) = b.K

"""
    FerriteProblems.getunknowns(p::FerriteProblem)

Get the current vector of unknowns from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FerriteProblems.getunknowns`
"""
FESolvers.getunknowns(b::FEBuffer) = b.x

"""
    FerriteProblems.getresidual(p::FerriteProblem)

Get the current residual vector from `p`. 
Note that this function belongs to `FESolvers.jl`,
but can be accessed via `FerriteProblems.getresidual`
"""
FESolvers.getresidual(b::FEBuffer) = b.r

"""
    FerriteProblems.getneumannforce(p::FerriteProblem)

Get the current external force vector caused by 
Neumann boundary conditions. Note that this vector 
does not include external forces added during the 
cell assembly; only forces added with the `NeumannHandler`
"""
get_external_force(b::FEBuffer) = b.f 

"""
    FerriteProblems.get_tolerance_scaling(p::FerriteProblem)

Get the `TolScaling` type that controls the convergence 
measure to be compared with the solver's tolerance
"""
get_tolerance_scaling(b::FEBuffer) = b.tolscaling


# Variables that will also be updated via special functions
# Unknowns 
"""
    FerriteProblems.get_old_unknowns(p::FerriteProblem)

Get the vector of unknowns from the previously converged step
"""
get_old_unknowns(b::FEBuffer) = b.xold 

"""
    FerriteProblems.update_unknowns!(p::FerriteProblem)

Update the vector of "old" unknowns to the values of the current vector of unknowns
"""
update_unknowns!(b::FEBuffer) = copy!(b.xold, b.x) # internal

# State variables
"""
    FerriteProblems.get_state(p::FerriteProblem)

Get the current state variables
"""
get_state(b::FEBuffer) = get_state(getassemblybuffer(b))

"""
    FerriteProblems.getoldstate(p::FerriteProblem)

Get the state variables from the previously converged step
"""
get_old_state(b::FEBuffer) = get_old_state(getassemblybuffer(b))

# Time
"""
    FerriteProblems.get_time(p::FerriteProblem)

Get the current time
"""
get_time(b::FEBuffer) = b.time

"""
    FerriteProblems.getoldtime(p::FerriteProblem)

Get time of the previous converged step
"""
get_old_time(b::FEBuffer) = b.old_time

"""
    FerriteProblems.settime!(p::FerriteProblem, new_time)

Set the current time to `new_time`
Called when starting a new step (or when attempting the same 
step number with a new time increment)
"""
set_time!(b::FEBuffer, new_time) = (b.time = new_time) # internal

"""
    FerriteProblems.update_time!(p::FerriteProblem)

Update the old time to the current time. 
Called after convergence
"""
update_time!(b::FEBuffer) = (b.old_time = get_time(b)) # internal
