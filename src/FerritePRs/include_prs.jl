# In this file useful PR's to Ferrite.jl can be added as a temporary solution
# Once merged, they should be included via a `@static` construct if the Ferrite version is old.

include("PR457.jl") # global_dof_range
include("PR427.jl") # add!(::ConstraintHandler, ::Dirichlet) for nonconcrete celltypes

@static if !isdefined(Ferrite, :apply_analytical!)
    include("PR532.jl")
end