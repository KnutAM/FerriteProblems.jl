"""
    FEDefinition(dh, ch, nh, cv, m, bl, initialstate, ic)
    FEDefinition(;dh, ch, cv, m,
        nh=NeumannHandler(dh), bl=nothing,
        initialstate=create_empty_states(dh,m),
        ic=(),
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
* `ic`: Initial conditions. NamedTuple with a function `f(x)` for each field that 
  has a nonzero initial condition. Used by the [`initial_conditions!`](@ref) function.
* `cc::ConvergenceCriterion`: Determines how to calculate the convergence measure including scaling

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
**Note**: *The following setup is not needed if you have one material*
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
struct FEDefinition{DH,CH,NH,CV,M,BL,ST,IC,CC}
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
    cc::CC   # Convergence criterion
end
function FEDefinition(;
    dh, ch, cv, m,          # Required for all simulations
    nh=NeumannHandler(dh),  # Can be just an empty handler
    bl=nothing,
    initialstate=create_empty_states(dh,m),
    ic=(),
    cc=AbsoluteResidual()
    )
    return FEDefinition(dh, ch, nh, cv, m, bl, initialstate, ic, cc)
end

# Note: The following functions are also overloaded for the entire ::FerriteProblem,
# and only this version is documented. 
"""
    FP.getdh(p::FerriteProblem)

Get `dh::Ferrite.AbstractDofHandler` from the `FEDefinition`
"""
getdh(def::FEDefinition) = def.dh

"""
    FP.getch(p::FerriteProblem)

Get the `ConstraintHandler` from the `FEDefinition`
"""
getch(def::FEDefinition) = def.ch

"""
    FP.getnh(p::FerriteProblem)

Get the `NeumannHandler` from the `FEDefinition`
"""
getnh(def::FEDefinition) = def.nh

"""
    FP.getcv(p::FerriteProblem)

Get the cell values from the `FEDefinition`. 
Note that this could also be a `Tuple` or `NamedTuple` depending on 
what was initially given to `FEDefinition`
"""
getcv(def::FEDefinition) = def.cv

"""
    FP.getmaterial(p::FerriteProblem)

Get the material from the `FEDefinition`
"""
getmaterial(def::FEDefinition) = def.m

"""
    FP.getbodyload(p::FerriteProblem)

Get the bodyload given to the `FEDefinition`
"""
getbodyload(def::FEDefinition) = def.bl


# Material Cache, default to nothing
allocate_material_cache(def::FEDefinition) = allocate_material_cache(getmaterial(def))
allocate_material_cache(materials::Dict) = Dict(key=allocate_material_cache(material) for (key,material) in materials)
allocate_material_cache(materials::Union{Tuple,NamedTuple}) = map(allocate_material_cache, materials)
"""
    FP.allocate_material_cache(material)

In case the material requires a cache to be available during the element routine,
this function can be overloaded for the specific material to define such a cache
to be included in the `FerriteAssembly.CellBuffer`
"""
allocate_material_cache(::Any) = nothing

# Default state creation (when not used)
"""
    FP.create_empty_states(dh, material)

Used to create empty states in case state variables aren't used in the simulation
by calling `FerriteAssembly.create_states` without any input (returning nothing for each cell)
"""
create_empty_states(dh::Ferrite.AbstractDofHandler, ::Any) = create_states(dh)  # internal
function create_empty_states(dh::Ferrite.AbstractDofHandler, m::Dict)
    try     
        return Dict(key=>create_states(dh, getcellset(dh, key)) for key in keys(m))
    catch e
        isa(e, KeyError) && println("If the material is a Dict, each key must correspond to a cellset in the grid")
        rethrow(e)
    end
end