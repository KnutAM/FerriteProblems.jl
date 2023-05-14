"""
Create the `FEDefinition` which can later be used to create a complete 
`FerriteProblem`. 

    FEDefinition(;dh, ch, material, cellvalues, 
        nh=NeumannHandler(dh),
        convergence_criterion=AbsoluteResidual(), 
        initial_conditions=NamedTuple(),
        user_data=nothing, 
        autodiffbuffer=Val(false), threading=Val(false)
        )

**Single domain:** This first definition sets up a single domain simulation 
(same `material`, `cellvalues`, interpolations, etc. everywhere). 

    FEDefinition(domains::Vector{<:AssemblyDomain}}; ch, 
        nh=NeumannHandler(dh),
        convergence_criterion=AbsoluteResidual(), 
        initial_conditions=NamedTuple(), 
        autodiffbuffer=Val(false), threading=Val(false)
        )
    
**Multiple domains:** This second definition setups up multi-domain simulations 
by using `FerriteAssembly.jl`'s `AssemblyDomain`. 
The `dh` used for the default NeumannHandler is extracted from `domains`. 

## Common inputs
For both constructors cases, the following inputs can be supplied
* `ch`: The constraint handler, `ConstraintHandler` (`Ferrite.jl`)
* `nh`: The neumann bc handler, `NeumannHandler` (`FerriteNeumann.jl`)
* `initial_conditions`: NamedTuple with a function `f(x)` for each field that has a nonzero 
  initial condition. Used by the `Ferrite.jl`'s `apply_analytical!` function.
  Example: `initial_conditions = (u = x -> Vec((x[1]/10, 0.0)), p = x -> -100*x[2])`. 
  For fields not given here, the initial condition is zeros everywhere.  
* `convergence_criterion`: Determines how to calculate the convergence measure including scaling.
  See [`ConvergenceCriterion`](@ref)
* `autodiffbuffer`: Should `FerriteAssembly.jl`'s `AutoDiffCellBuffer` be used? 
  This will make the assembly faster if automatic differentiation is used, and can also be used 
  without automatic differentiation (but with a slight extra computational overhead)
  `Bool` input can also be given (will be converted internally to `Val`)
* `threading`: Should threading be used? 
  `Bool` input can also be given (will be converted internally to `Val`)
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