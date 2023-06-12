"""
    CustomStiffness(material, stiffness_type::Symbol)

`CustomStiffness` is a wrapper for a `material`, which makes it possible to define custom behavior 
without changing the type of the material. This wrapper is intended to change how the 
element stiffness, `Ke`, is calculated, by branching on the `stiffness_type` value. 
This value can be changed by calling `set_stiffness_type!(::CustomStiffness, ::Symbol)`. 
To use this wrapper, overload the `element_routine!` for `::CustomStiffness{<:MyMat}` as follows. 
```julia
import FerriteAssembly as FA
import FerriteProblems as FP
function FA.element_routine!(Ke, re, new_state, ae, m::CustomStiffness{<:MyMat}, args...; kwargs...)
    if FP.get_stiffness_type(m) == :true 
        FA.element_routine!(Ke, re, new_state, ae, m.material, args...; kwargs...)
    elseif FP.get_stiffness_type(m) == :modified_picard # Call an alternative element routine
        FA.element_routine!(Ke, re, new_state, ae, ModifiedPicardWrapper(m.material), args...; kwargs...)
    elseif FP.get_stiffness_type(m) == :relaxed
        FA.element_routine!(Ke, re, new_state, ae, m.material, args...; kwargs...)
        relax_stiffness!(Ke, m.material) # User routine to modify the stiffness
    else
        throw(ArgumentError(join(("jacobian_type: ", get_stiffness_type(m), " is not supported"))))
    end
end
```
where `ModifiedPicardWrapper` is a wrapper to change how the stiffness is calculated, i.e.
`FA.element_routine!(Ke, re, new_state, ae, m::ModifiedPicardWrapper{<:MyMat}, args...; kwargs...)`
is implementated. The `CustomStiffness` wrapper uses the `FerriteAssembly.unwrap_material_for_ad` 
feature to support calculating the stiffness for `m.material` using automatic differentiation as usual.
It also supports the `set_jacobian_type` method, such that the type requested by `FESolvers.UpdateSpec`
can be respected by the `CustomStiffness`-wrapped material. 
"""
mutable struct CustomStiffness{M}
    const material::M
    stiffness_type::Symbol
end
function FerriteAssembly.unwrap_material_for_ad(cs::CustomStiffness)
    return FerriteAssembly.unwrap_material_for_ad(cs.material)
end
get_stiffness_type(cs::CustomStiffness) = cs.stiffness_type

function FerriteAssembly.element_routine!(Ke, re, new_state, ae, m::CustomStiffness{M}, args...; kwargs...) where M
    msg = join(("You must implement element_routine!, normally for m::CustomStiffness{<:", nameof(M), "}"))
    throw(ArgumentError(msg))
end

function FerriteAssembly.element_residual!(re, new_state, ae, m::CustomStiffness, args...; kwargs...)
    FerriteAssembly.element_residual!(re, new_state, ae, m.material, args...; kwargs...)
end

function FerriteAssembly.create_cell_state(m::CustomStiffness, args...; kwargs...)
    FerriteAssembly.create_cell_state(m.material, args...; kwargs...)
end

function FerriteAssembly.allocate_cell_cache(m::CustomStiffness, args...)
    FerriteAssembly.allocate_cell_cache(m.material, args...)
end

# Overload `set_jacobian_type` for `CustomStiffness`:
function set_jacobian_type(m::CustomStiffness, type::Symbol)
    m.stiffness_type = type
    return m
end

# Do nothing of no type is specified (standard solvers):
set_jacobian_type(m::CustomStiffness, ::Nothing) = m 

# Throw error if non-supported jacobian type is given:
function set_jacobian_type(::CustomStiffness, _)
    throw(ArgumentError("Only type::Symbol is supported for CustomStiffness"))
end