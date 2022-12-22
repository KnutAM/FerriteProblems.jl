# https://github.com/Ferrite-FEM/Ferrite.jl/pull/427
function Ferrite.add!(ch::ConstraintHandler{<:MixedDofHandler}, dbc::Dirichlet)
    dbc_added = false
    for fh in ch.dh.fieldhandlers
        if _overlaps(fh, dbc)
            # If dbc have dofs not in fh, then these will be removed from dbc, hence deepcopy
            # In this case, add! will warn, unless warn_check_cellset=false
            #add!(ch, fh, deepcopy(dbc), warn_check_cellset=false)
            add!(ch, fh, deepcopy(dbc)) # supressing warning not easy except from inside pr
            dbc_added = true
        end
    end
    dbc_added || @warn("No overlap between dbc::Dirichlet and fields in the ConstraintHandler's MixedDofHandler")
    return ch
end

function _overlaps(fh::FieldHandler, dbc::Dirichlet)
    dbc.field_name in Ferrite.getfieldnames(fh) || return false # Must contain the correct field
    for (cellid, _) in dbc.faces
        cellid in fh.cellset && return true
    end
    return false
end