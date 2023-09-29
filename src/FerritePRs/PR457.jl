# Based on #457: Add global_dof_range to get the global index range for a field.
function global_dof_range(dh::DofHandler, field_name::Symbol)
    dofs = Set{Int}()
    for sdh in dh.subdofhandlers
        if field_name ∈ Ferrite.getfieldnames(sdh)
            _global_dof_range!(dofs, sdh, field_name)
        end
    end
    return sort!(collect(Int, dofs))
end
function _global_dof_range!(dofs, sdh::SubDofHandler, field_name)
    cellset = getcellset(sdh)
    eldofs = celldofs(dh, first(cellset))
    field_range = dof_range(sdh, field_name)
    for i in cellset
        celldofs!(eldofs, dh, i)
        for j in field_range
            @inbounds d = eldofs[j]
            d in dofs || push!(dofs, d)
        end
    end
end
 