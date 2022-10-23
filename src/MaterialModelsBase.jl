# Mechanical materials that are defined according to the `MaterialModelsBase.jl`
# have already built-in support for convenience
get_material_cache(m::AbstractMaterial) = get_cache(m)

getstatefun(m::AbstractMaterial) = _->initial_material_state(m)

function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, state::Vector{<:AbstractMaterialState},
    ae::AbstractVector, material::AbstractMaterial, cellvalues::CellVectorValues, 
    dh_fh, Δt, buffer::CellBuffer)
    cache = buffer.cache
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        # For each integration point, compute stress and material stiffness
        ϵ = function_symmetric_gradient(cellvalues, q_point, ae) # Total strain
        σ, D, state[q_point] = material_response(material, ϵ, state[q_point], Δt, cache)

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            ∇δN = shape_symmetric_gradient(cellvalues, q_point, i)
            re[i] += (∇δN ⊡ σ) * dΩ # add internal force to residual
            ∇δN_D = ∇δN ⊡ D         # temporary value for speed
            for j in 1:n_basefuncs
                ∇N = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += (∇δN_D ⊡ ∇N) * dΩ
            end
        end
    end
end;
