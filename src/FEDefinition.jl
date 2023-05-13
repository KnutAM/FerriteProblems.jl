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
struct FEDefinition{DH,CH,NH,CC}
    # FE-handlers
    dh::DH  # DofHandler/MixedDofHandler
    ch::CH  # ConstraintHandler
    nh::NH  # NeumannHandler
    # Convergence criterion
    convergence_criterion::CC
    # For construction of FEBuffer, should not be used for performance critical code
    initial_conditions::NamedTuple  # NamedTuple: (fieldname=f(x),)
    autodiffbuffer::Val
    threading::Val
    domains::Union{<:AssemblyDomain, <:Vector{<:AssemblyDomain}}
end
_extract_dofhandler(domain::AssemblyDomain) = domain.sdh.dh
_extract_dofhandler(domains::Vector{<:AssemblyDomain}) = domains[1].sdh.dh
function FEDefinition(domain::Union{AssemblyDomain, Vector{<:AssemblyDomain}}; 
        ch, nh=NeumannHandler(_extract_dofhandler(domain)),
        convergence_criterion=AbsoluteResidual(), initial_conditions=NamedTuple(), 
        autodiffbuffer=Val(false), threading=Val(false)
        )
    dh = _extract_dofhandler(domain)
    return FEDefinition(dh, ch, nh, convergence_criterion, initial_conditions, autodiffbuffer, threading, domain)
end
# Constructor without creating AssemblyDomain first, used for single domains
function FEDefinition(;
    dh, ch, material, cellvalues, # Required for all simulations
    user_data=nothing, kwargs...)
    domain = AssemblyDomain("full", dh, material, cellvalues; user_data=user_data)
    return FEDefinition(domain; ch=ch, kwargs...)
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