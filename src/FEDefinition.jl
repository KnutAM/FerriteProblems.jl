"""
    FEDefinition(dh, ch, nh, cv, m, bl, initialstate, ic)
    FEDefinition(;
        dh, ch, cv, m,
        nh=NeumannHandler(dh), bl=nothing, ic=(), 
        initialstate=_create_states(dh,m,cv,ic),
        cc=AbsoluteResidual(), colors=nothing,
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
* `ic`: Initial conditions. NamedTuple with a function `f(x)` for each field that has a nonzero 
  initial condition. Used by the `Ferrite.jl`'s `apply_analytical!` function 
  (included here if Ferrite version too old).
* `initialstate`: The initial state variables for each cell in the grid. By default, this is 
  automatically created by `FerriteAssembly.create_states` with customization provided by
  overloading `FerriteAssembly.create_cell_state` and is not necessary to supply here. 
* `cc::ConvergenceCriterion`: Determines how to calculate the convergence measure including scaling
* `colors`: Supply to use threaded assembly, can be created by `Ferrite.create_coloring`

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
struct FEDefinition{DH,CH,NH,CV,M,BL,ST,IC,CC,COLOR}
    # FE-handlers
    dh::DH  # DofHandler/MixedDofHandler
    ch::CH  # ConstraintHandler
    nh::NH  # NeumannHandler
    # Items related to each cell type in the grid
    # If multiple cell types, these should be tuples
    cv::CV  # cellvalues
    m::M    # material
    bl::BL  # body loads/source terms
    # Initial values
    ic::IC  # NamedTuple: (fieldname=f(x),)
    initialstate::ST
    # Convergence criterion
    cc::CC   
    # Threaded assembly if colors!=nothing
    colors::COLOR
end
function FEDefinition(;
    dh, ch, cv, m,          # Required for all simulations
    nh=NeumannHandler(dh),  # Can be just an empty handler
    bl=nothing,
    ic=(),
    initialstate=_create_states(dh,m,cv,ic),
    cc=AbsoluteResidual(),
    colors=nothing
    )
    return FEDefinition(dh, ch, nh, cv, m, bl, ic, initialstate, cc, colors)
end

# Note: The following functions are also overloaded for the entire ::FerriteProblem,
# and only this version is documented. 
"""
    FerriteProblems.getdh(p::FerriteProblem)

Get `dh::Ferrite.AbstractDofHandler` from the `FerriteProblem`
"""
getdh(def::FEDefinition) = def.dh

"""
    FerriteProblems.getch(p::FerriteProblem)

Get the `ConstraintHandler` from the `FerriteProblem`
"""
getch(def::FEDefinition) = def.ch

"""
    FerriteProblems.getnh(p::FerriteProblem)

Get the `NeumannHandler` from the `FerriteProblem`
"""
getnh(def::FEDefinition) = def.nh

"""
    FerriteProblems.getcv(p::FerriteProblem)

Get the cell values from the `FerriteProblem`. 
Note that this could also be a `Tuple` or `NamedTuple` depending on 
what was initially given to `FEDefinition`
"""
getcv(def::FEDefinition) = def.cv

"""
    FerriteProblems.getmaterial(p::FerriteProblem)

Get the material from the `FerriteProblem`
"""
getmaterial(def::FEDefinition) = def.m

"""
    FerriteProblems.getbodyload(p::FerriteProblem)

Get the bodyload given to the `FerriteProblem`
"""
getbodyload(def::FEDefinition) = def.bl


# Material Cache, default to nothing
"""
    FerriteProblems.allocate_material_cache(material, cellvalues)

In case the material requires a cache to be available during the element routine,
this function can be overloaded for the specific material to define such a cache
to be included in the `FerriteAssembly.CellBuffer`
"""
allocate_material_cache(args...) = nothing

# Top level call from definition
allocate_material_cache(def::FEDefinition) = allocate_material_cache(getmaterial(def), getcv(def))

# Multiple materials
function allocate_material_cache(materials::Dict, cellvalues)
    mtrl_keys = keys(materials)
    cv_ = FerriteAssembly._makedict(cellvalues, mtrl_keys)
    return Dict(key=allocate_material_cache(material[key], cv_[key]) for key in mtrl_keys)
end

# MixedDofHandler
function allocate_material_cache(materials::Tuple, cellvalues)
    cv_ = FerriteAssembly._maketuple(cellvalues, length(materials))
    return map(allocate_material_cache, materials, cv_)
end

"""
    _create_states(dh::AbstractDofHandler, material, cellvalues, initial_conditions)

Create the state variables by calling `FerriteAssembly.create_states` after calculating 
the degree of freedom values with [`FerriteProblems.initial_conditions!`](@ref) and the 
`initial_conditions` input.  
"""
function _create_states(dh, material, cellvalues, initial_conditions)
    a = zeros(ndofs(dh))
    if length(initial_conditions) > 0
        foreach(ic->initial_conditions!(a, dh, ic[1], ic[2]), pairs(initial_conditions))
    end
    create_states(dh, material, cellvalues, a)
end

"""
    dothreaded(def::FEDefinition)

Trait-based dispatch to determine if threaded assembly should be used.
Returns a `Val{t::Bool}`. 
"""
dothreaded(def::FEDefinition) = _dothreaded(def.colors)
_dothreaded(::Nothing) = Val{false}()
_dothreaded(::Any) = Val{true}()
