"""
    FEDefinition(
        domainspec::Union{DomainSpec,Dict{String,DomainSpec}};
        ch, lh=LoadHandler(dh),
        convergence_criterion=AbsoluteResidual(), 
        initial_conditions=NamedTuple(),
        )

Create the `FEDefinition` which can later be used to 
create a complete `FerriteProblem`. `domainspec` is
the `FerriteAssembly` domain specification.
 
In addition, the following keyword arguments can be given
* `ch`: The constraint handler, `ConstraintHandler` (`Ferrite.jl`)
* `lh`: The external load handler, `LoadHandler` (`FerriteAssembly.jl`)
* `initial_conditions`: NamedTuple with a function `f(x)` for each field that has a nonzero 
  initial condition. Used by the `Ferrite.jl`'s `apply_analytical!` function.
  Example: `initial_conditions = (u = x -> Vec((x[1]/10, 0.0)), p = x -> -100*x[2])`. 
  For fields not given here, the initial condition is zeros everywhere.  
* `convergence_criterion`: Determines how to calculate the convergence measure including scaling.
  See [`ConvergenceCriterion`](@ref)
"""
struct FEDefinition{DH,CH,LH,CC}
    # FE-handlers
    dh::DH  # DofHandler
    ch::CH  # ConstraintHandler
    lh::LH  # LoadHandler
    # Convergence criterion
    convergence_criterion::CC
    # For construction of FEBuffer, should not be used for performance critical code
    initial_conditions::NamedTuple  # NamedTuple: (fieldname=f(x),)
    domains # DomainSpec or Dict{String, DomainSpec}
end
_extract_dofhandler(domain::DomainSpec) = domain.sdh.dh
_extract_dofhandler(domains::Dict{String, DomainSpec}) = first(values(domains)).sdh.dh
function FEDefinition(domain::Union{DomainSpec, Dict{String, DomainSpec}}; 
        ch, lh=LoadHandler(_extract_dofhandler(domain)),
        convergence_criterion=AbsoluteResidual(), initial_conditions=NamedTuple(), 
        )
    dh = _extract_dofhandler(domain)
    return FEDefinition(dh, ch, lh, 
        convergence_criterion, 
        initial_conditions, 
        domain)
end

# Note: The following functions are also overloaded for the entire ::FerriteProblem,
# and only this version is documented. 
"""
    FerriteProblems.get_dofhandler(p::FerriteProblem)

Get `dh::Ferrite.AbstractDofHandler` from the `FerriteProblem`
(Technically overloaded from `FerriteAssembly`, but accessible via `FerriteProblems`)
"""
get_dofhandler(def::FEDefinition) = def.dh

"""
    FerriteProblems.get_constrainthandler(p::FerriteProblem)

Get the `ConstraintHandler` from the `FerriteProblem`
"""
get_constrainthandler(def::FEDefinition) = def.ch

"""
    FerriteProblems.get_loadhandler(p::FerriteProblem)

Get the `LoadHandler` from the `FerriteProblem`
"""
get_loadhandler(def::FEDefinition) = def.lh

