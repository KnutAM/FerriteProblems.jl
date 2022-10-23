struct FEDefinition{DH,CH,NH,CV,M,BL,ST,IC}
    # FE-handlers
    dh::DH  # DofHandler/MixedDofHandler
    ch::CH  # ConstraintHandler
    nh::NH  # NeumannHandler
    # Items related to each cell type in the grid
    # If multiple cell types, these should be tuples
    cv::CV  # cellvalues
    m::M    # material
    bl::BL  # body loads/source terms
    initialstate::ST
    # Initial value of dofs 
    ic::IC  # NamedTuple: (fieldname=f(x),)
end
function FEDefinition(;
    dh, ch, cv, m,          # Required for all simulations
    nh=NeumannHandler(dh),  # Can be just an empty handler
    bl=nothing,
    initialstate=create_empty_states(dh,m),
    ic=()
    )
    return FEDefinition(dh, ch, nh, cv, m, bl, initialstate, ic)
end
"""
    FEDefinition(dh, ch, nh, cv, m, bl, initialstate)
    FEDefinition(;dh, ch, cv, m,
        nh=NeumannHandler(dh), bl=nothing,
        initialstate=create_empty_states(dh,m)
        )

The definition of the FE-problem, this data describes the finite 
element problem as a function of time. For how long to run and how 
to step through time is typically done by `FESolvers.jl`. 

## Fields/inputs
See more information below for items marked with *

* `dh`: The dof handler, `Union{DofHandler,MixedDofHandler}` (`Ferrite.jl`)
* `ch`: The constraint handler, `ConstraintHandler` (`Ferrite.jl`)
* `nh`: The neumann bc handler, `NeumannHandler` (`FerriteNeumann.jl`)
* `cv`: The cellvalues, e.g. `CellValues` (*)
* `m`: The material definition, user defined type - passed into element. (*)
* `bl`: Source term/body load, user defined type - available from the element routine (*)
* `initialstate`: The initial state variables for each cell in the grid, 
  can created by `FerriteAssembly.create_states` (*)

## `MixedDofHandler`
When the `MixedDofHandler` is used, we have an outer loop over each of its `FieldHandler`s.
Therefore, `cv`, `m`, and `bl` may be provided for each of these fields by providing them as 
a tuple. Otherwise, they are duplicated for each field (by reference)

`initialstate` must be a tuple containing `Dict{Int}` with states for the cells in the relevant 
`FieldHandler`'s cellset. 

**Multiple `CellValue`s for each element**
In coupled problems, each element might require multiple `CellValue`s. 
If the same should be used for all fields (even if only one `FieldHandler`), 
a `NamedTuple` is the easiest (the alternative is a `Tuple` of `Tuple`s). 
If different multiple `CellValues` should be used for each field, then pass as 
`Tuple` of `NamedTuple`/`Tuple`

## Multiple materials
***Note***: *The following setup is not needed if you have one material*
*for each `FieldHandler` when using a `MixedDofHandler`*

If you want to have different materials on different parts of the grid, 
you should first define those cellsets, and then define each material 
as a value of `m::Dict{String}` with the name of the cellset as key. 

You may do the same for `bl` and `cv`. 
Otherwise, the same values will be used for every cellset. 

**Combine with `MixedDofHandler`**: 
Follow the required input for `MixedDofHandler`, 
but wrap items (at least the material and `initialstate`)
in `Dict`s as described above

**State variables**: 
`initialstate::Dict{String}` must mirror the datastructure of `m::Dict{String}`:
For each `key` in `m::Dict{String}`, `initialstate[key]` holds a collection of states 
for the cells in `getcellset(dh,key)`.

If `MixedDofHandler`, then `initialstate[key]::NTuple{N,<:Dict{Int}}`, 
and if `DofHandler`, then `initialstate[key]::Dict{Int}`
"""

getdh(def::FEDefinition) = def.dh
getch(def::FEDefinition) = def.ch
getnh(def::FEDefinition) = def.nh
getcv(def::FEDefinition) = def.cv
getmaterial(def::FEDefinition) = def.m
getbodyload(def::FEDefinition) = def.bl
getic(def::FEDefinition) = def.ic

# Material Cache, default to nothing
get_material_cache(def::FEDefinition) = get_material_cache(getmaterial(def))
get_material_cache(materials::Dict) = Dict(key=get_material_cache(material) for (key,material) in materials)
get_material_cache(materials::Union{Tuple,NamedTuple}) = map(get_material_cache, materials)
get_material_cache(::Any) = nothing

# Default state creation (when not used)
create_empty_states(dh::Ferrite.AbstractDofHandler, ::Any) = create_states(dh)
function create_empty_states(dh::Ferrite.AbstractDofHandler, m::Dict)
    try     
        return Dict(key=>create_states(dh, getcellset(dh, key)) for key in keys(m))
    catch e
        isa(e, KeyError) && println("If the material is a Dict, each key must correspond to a cellset in the grid")
        rethrow(e)
    end
end