# 
# #457: Add global_dof_range to get the global index range for a field.
function global_dof_range(dh::MixedDofHandler, field_name::Symbol)
    dofs = Set{Int}()
    for fh in dh.fieldhandlers
        if field_name ∈ Ferrite.getfieldnames(fh)
            _global_dof_range!(dofs, dh, fh, field_name, fh.cellset)
        end
    end
    return sort!(collect(Int, dofs))
end
function global_dof_range(dh::DofHandler, field_name::Symbol)
    dofs = Set{Int}()
    _global_dof_range!(dofs, dh, dh, field_name, 1:getncells(dh.grid))
    return sort!(collect(Int, dofs))
end
function _global_dof_range!(dofs, dh, dh_fh, field_name, cellset)
    eldofs = celldofs(dh, first(cellset))
    field_range = dof_range(dh_fh, field_name)
    for i in cellset
        celldofs!(eldofs, dh, i)
        for j in field_range
            @inbounds d = eldofs[j]
            d in dofs || push!(dofs, d)
        end
    end
end
#
#
# 