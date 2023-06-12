"""
    FEBuffer(K,x,r,f,xold,cellbuffer,state,old_state,time,old_time)

A buffer to hold all values that are required to simulate, 
but that are uniqely defined from the simulation definition
"""
mutable struct FEBuffer{T,KT<:AbstractMatrix{T},AB,ST,TS} # internal 
    const K::KT
    const x::Vector{T}
    const r::Vector{T}
    const f::Vector{T}
    const xold::Vector{T}
    const assembly_buffer::AB   # Buffer from FerriteAssembly.setup_assembly
    const state::ST             # Also following output 
    const old_state::ST         # from setup_assembly
    time::T
    old_time::T
    const tolscaling::TS
end

function FEBuffer(def::FEDefinition)
    n = ndofs(getdh(def))
    K = create_sparsity_pattern(getdh(def), getch(def))
    x, r, f = zeros.((n,n,n))
    ic = def.initial_conditions # ::NamedTuple 
    foreach((field_name, fun)->apply_analytical!(x, getdh(def), field_name, fun), keys(ic), values(ic))
    xold = deepcopy(x)
    buffer, new_state, old_state = _setup_assembly(def, def.domains; a=x)
    time = 0.0
    old_time = 0.0
    tol_scaling = TolScaling(def.convergence_criterion, def)
    return FEBuffer(K, x, r, f, xold, buffer, new_state, old_state, time, old_time, tol_scaling)
end

function _setup_assembly(def::FEDefinition, d::AssemblyDomain; a)
    return setup_assembly(d.sdh, d.material, d.cellvalues; 
        a=a, threading=def.threading, autodiffbuffer=def.autodiffbuffer,
        cellset=d.cellset, colors=d.colors, user_data=d.user_data
        )
end
function _setup_assembly(def::FEDefinition, d::Vector{<:AssemblyDomain}; a)
    return setup_assembly(d; a=a, threading=def.threading, autodiffbuffer=def.autodiffbuffer)
end

# Standard get functions
getassemblybuffer(b::FEBuffer) = b.assembly_buffer
"""
    FerriteProblems.get_material(p::FerriteProblem)
    FerriteProblems.get_material(p::FerriteProblem, domain_name::String)

Get the material in `p`. For multiple domains, it is necessary to give the `domain_name`
for where to get the material. Note that this is type-unstable and should be avoided in 
performance-critical code sections. This function belongs to `FerriteAssembly.jl`, 
but can be accessed via `FerriteProblems.get_material`.
"""
FerriteAssembly.get_material(b::FEBuffer{<:Any,<:Any,<:FerriteAssembly.AbstractDomainBuffer}) = get_material(getassemblybuffer(b))
FerriteAssembly.get_material(b::FEBuffer{<:Any,<:Any,<:Dict{String}}, name::String) = get_material(getassemblybuffer(b), name)
function FerriteAssembly.get_material(::FEBuffer{<:Any,<:Any,<:Dict{String}})
    throw(ArgumentError("get_material requires the domain name for multiple domain simulations"))
end

"""
    set_jacobian_type(material, type)

Define this function for your material, and the type given by `FESolvers.UpdateSpec` 
when the solver asks for different jacobian calculation. The function should return the 
modified material according to the given type. See `CustomStiffness` for an example, or 
use this wrapper directly if that is sufficient. 
"""
set_jacobian_type(material, type) = material # Default to not implemented, silently doing nothing if not supported.

function change_material_jacobian_type!(b::FEBuffer, type)
    FerriteAssembly.modify_material!(m->set_jacobian_type(m, type), getassemblybuffer(b))
end

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
getneumannforce(b::FEBuffer) = b.f 

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
