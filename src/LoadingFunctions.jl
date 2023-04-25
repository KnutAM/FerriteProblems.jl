# To be able to reliably save the definition, functions passed to `ConstraintHandler`, `NeumannHandler`, 
# and for initial conditions must be defined as actual functions that are also available upon reloading. 
# In this file, a few common functions can be defined to simplify this.
# In general, use `Returns(value)` for all constant values

# Scalar functions
struct RampHold{T} <: Function
    t_rise::T
    val_change::T
    val_initial::T
end
RampHold(;t_rise, val_change, val_initial) = RampHold(t_rise, val_change, val_initial)

(rh::RampHold{T})(x, t, args...) where {T} = rh.val_initial + rh.val_change*clamp(t/rh.t_rise, 0, 1)

# Vector functions
struct VectorFunction{Dim,FT<:Tuple} <: Function
    funs::FT
end
VectorFunction(args...) = VectorFunction(args)
VectorFunction(arg) = VectorFunction((arg,))
function VectorFunction(args::Tuple)
    _make_function(arg) = arg
    _make_function(arg::Number) = Returns(arg)
    funs = map(_make_function, args)
    return VectorFunction{length(funs), typeof(funs)}(funs)
end
(vf::VectorFunction{1})(args...) = Vec(vf.funs[1](args...))
(vf::VectorFunction{2})(args...) = Vec(vf.funs[1](args...), vf.funs[2](args...))
(vf::VectorFunction{3})(args...) = Vec(vf.funs[1](args...), vf.funs[2](args...), vf.funs[3](args...))
