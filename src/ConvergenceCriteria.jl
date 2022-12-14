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
* `assemscaling`: Used in `FerriteAssembly.doassemble!` to give potential scaling contribution based on each element's residual
* `buffer`: Used to pre-calculate values, such as dof-ranges for each field that is used when calculating the convergence measure. 
"""
struct TolScaling{C,S,B}
    criterion::C    # Parameters and scaling type: assemscaling and buffer can be created from this one
    assemscaling::S # Scaling to be given to doassemble!
    buffer::B       # Buffer to hold precalculated values, such as dofs for each field
end
TolScaling(criterion, def) = TolScaling(criterion, make_assemscaling(criterion, def), nothing)

"""
    make_assemscaling(criterion, def, threaded::Val{false})

Create the actual assemscaling for the given `criterion`.

    make_assemscaling(criterion, def, threaded::Val{true})

If `threaded::Val{true}`, then call `make_assemscaling(criterion, def, ::Val{false})` 
for each thread to support threaded assembly. Note that the function has a definition 
`make_assemscaling(criterion,def)=make_assemscaling(criterion,def,dothreaded(def))` to simplify calling it.
If required, a `criterion` should overload `make_assemscaling(criterion,def,::Val{false})`.
The default, if not overloaded, returns `nothing`. 
That suffices if nothing needs to be calculated for each cell during assembly. 
"""
make_assemscaling(criterion, def) = make_assemscaling(criterion, def, dothreaded(def))
make_assemscaling(::Any, _, ::Val{false}) = nothing
make_assemscaling(criterion, def, ::Val{true}) = create_threaded_scalings(make_assemscaling(criterion, def, Val{false}()))

"""
    FESolvers.calculate_convergence_measure(::ConvergenceCriterion, scaling::TolScaling, r::Vector, ??a, iter, p::FerriteProblem)

Specialize this function for each `ConvergenceCriterion`. Note this function is called via the definition
```
calculate_convergence_measure(ts::TolScaling, args...) = calculate_convergence_measure(ts.criterion, ts, args...)
``` 
"""
FESolvers.calculate_convergence_measure(ts::TolScaling, args...) = FESolvers.calculate_convergence_measure(ts.criterion, ts, args...)

"""
    AbsoluteResidual()

The default convergence criterion that calculates the convergence measure as 
`???(sum([r[i]^2 for i ??? free dofs])` without any scaling. 
"""
struct AbsoluteResidual <: ConvergenceCriterion end

function FESolvers.calculate_convergence_measure(::AbsoluteResidual, ts, r, ??a, iter, p)
    sqrt(sum(i->r[i]^2, Ferrite.free_dofs(getch(p))))
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
return ???val
```
where the same `minfactor` is used for all fields if only a scalar value is given. 
"""
struct RelativeResidualElementScaling{T,F<:Union{AbstractFloat,Dict{Symbol},NamedTuple}} <: ConvergenceCriterion
    p::T
    minfactors::F
end 
RelativeResidualElementScaling(;p=2, minfactors=eps()) = RelativeResidualElementScaling(p, minfactors)

make_assemscaling(criterion::RelativeResidualElementScaling, def, ::Val{false}) = ElementResidualScaling(getdh(def), criterion.p)

function TolScaling(criterion::RelativeResidualElementScaling, def)
    dh = getdh(def)
    assemscaling = make_assemscaling(criterion, def)
    fdofs = Ferrite.free_dofs(getch(def))
    buffer = Dict(key=>intersect!(global_dof_range(dh, key), fdofs) for key in Ferrite.getfieldnames(dh))
    return TolScaling(criterion, assemscaling, buffer)
end

function FESolvers.calculate_convergence_measure(cc::RelativeResidualElementScaling, ts, r, args...)
    val = zero(eltype(r))
    factors = sum(ts.assemscaling).factors
    getminfactor(minfactors, args...) = minfactors
    getminfactor(minfactors::Union{Dict{Symbol},NamedTuple}, fieldname::Symbol) = minfactors[fieldname]
    minfactors = cc.minfactors
    for (key, dofs) in ts.buffer 
        factor = max(factors[key], getminfactor(minfactors, key))
        val += sum(i->(r[i]/factor)^2, dofs)    # Inside loop for better accuracy
    end
    return sqrt(val)
end
