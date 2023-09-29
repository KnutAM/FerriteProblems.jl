# Scaling of the tolerance for the finite element problem
"""
    ConvergenceCriterion

The abstract type `ConvergenceCriterion` is the supertype for all 
convergence criteria. 
"""
abstract type ConvergenceCriterion end

"""
    TolScaling(criterion::ConvergenceCriterion, def::FEDefinition)

The `TolScaling` type contains the `criterion` that determines how to scale the residuals to determine 
convergence. The constructor is specialized on `typeof(criterion)`, creating the following fields:
* `assemscaling`: `scaling` to be used by `FerriteAssembly` to give potential scaling contribution based on each element's residual
* `buffer`: Used to pre-calculate values, such as dof-ranges for each field that is used when calculating the convergence measure. 
"""
struct TolScaling{C,S,B}
    criterion::C    # Parameters and scaling type: assemscaling and buffer can be created from this one
    assemscaling::S # Scaling to be given to doassemble!
    buffer::B       # Buffer to hold precalculated values, such as dofs for each field
end
TolScaling(criterion, def) = TolScaling(criterion, make_assemscaling(criterion, def), nothing)

"""
    make_assemscaling(criterion, def)

Create the `scaling` for use in `FerriteAssembly` if required by the given `criterion`.
The default, if not overloaded, returns `nothing`. 
"""
make_assemscaling(criterion, def) = FerriteAssembly.NoScaling()

"""
    FESolvers.calculate_convergence_measure(::ConvergenceCriterion, scaling::TolScaling, r::Vector, Δa, iter, p::FerriteProblem)

Specialize this function for each `ConvergenceCriterion`. Note this function is called via the definition
```
calculate_convergence_measure(ts::TolScaling, args...) = calculate_convergence_measure(ts.criterion, ts, args...)
``` 
"""
FESolvers.calculate_convergence_measure(ts::TolScaling, args...) = FESolvers.calculate_convergence_measure(ts.criterion, ts, args...)

"""
    AbsoluteResidual()

The default convergence criterion that calculates the convergence measure as 
`√(sum([r[i]^2 for i ∈ free dofs])` without any scaling. 
"""
struct AbsoluteResidual <: ConvergenceCriterion end

function FESolvers.calculate_convergence_measure(::AbsoluteResidual, ts, r, Δa, iter, p)
    sqrt(sum(i->r[i]^2, Ferrite.free_dofs(get_constrainthandler(p))))
end

"""
    RelativeResidualElementScaling(p, minfactors::Union{AbstractFloat,NamedTuple}=eps())

Use `Ferriteassembly.ElementResidualScaling` with the exponent `p` to calculate the 
scaling for each field individually, based on the L2-norm of each cell's residual. 
To avoid issues when all cells have zero residual (e.g. in the first time step),
supply `minfactors` as the minimum scaling factor. 
The convergence measure is calculated with the following pseudo-code
```julia
val = 0.0
for field in Ferrite.getfieldnames(dh)
    factor = max(element_residual_scaling[field], minfactors[field])
    dofs = free_field_dofs[field]   # Get the non-constrained dofs for `field`
    val += sum(i->(r[i]/factor)^2, dofs)
end
return √val
```
where the same `minfactor` is used for all fields if only a scalar value is given. 
"""
struct RelativeResidualElementScaling{T,F<:Union{AbstractFloat,Dict{Symbol},NamedTuple}} <: ConvergenceCriterion
    p::T
    minfactors::F
end 
RelativeResidualElementScaling(;p=Val(2), minfactors=eps()) = RelativeResidualElementScaling(p, minfactors)

make_assemscaling(criterion::RelativeResidualElementScaling, def) = ElementResidualScaling(get_dofhandler(def), criterion.p)

function TolScaling(criterion::RelativeResidualElementScaling, def)
    dh = get_dofhandler(def)
    assemscaling = make_assemscaling(criterion, def)
    fdofs = Ferrite.free_dofs(get_constrainthandler(def))
    buffer = Dict(key=>intersect!(global_dof_range(dh, key), fdofs) for key in Ferrite.getfieldnames(dh))
    return TolScaling(criterion, assemscaling, buffer)
end

function FESolvers.calculate_convergence_measure(cc::RelativeResidualElementScaling, ts, r, args...)
    val = zero(eltype(r))
    factors = ts.assemscaling.factors
    getminfactor(minfactors, args...) = minfactors
    getminfactor(minfactors::Union{Dict{Symbol},NamedTuple}, fieldname::Symbol) = minfactors[fieldname]
    minfactors = cc.minfactors
    for (key, dofs) in ts.buffer 
        factor = max(factors[key], getminfactor(minfactors, key))
        val += sum(i->(r[i]/factor)^2, dofs)    # Inside loop for better accuracy
    end
    return sqrt(val)
end
