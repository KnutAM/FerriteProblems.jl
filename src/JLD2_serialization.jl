struct SerializedSubDofHandler
    cellset::Set{Int}
    fields::Vector{<:Pair{Symbol,<:Interpolation}}
end

struct SerializedDofHandler
    grid
    subdofhandlers::Vector{SerializedSubDofHandler}
end

JLD2.writeas(::Type{<:Ferrite.DofHandler}) = SerializedDofHandler
function JLD2.convert(::Type{<:Ferrite.DofHandler}, serialized_dh::SerializedDofHandler)
    dh = Ferrite.DofHandler(serialized_dh.grid)
    for serialized_sdh in serialized_dh.subdofhandlers
        sdh = Ferrite.SubDofHandler(dh, serialized_sdh.cellset)
        for (field, ip) in serialized_sdh.fields
            add!(sdh, field, ip)
        end
    end
    close!(dh)
    return dh
end
function JLD2.convert(::Type{SerializedDofHandler}, dh::Ferrite.DofHandler)
    subdofhandlers = SerializedSubDofHandler[]
    for sdh in dh.subdofhandlers
        fields = Pair{Symbol,Interpolation}[]
        for name in Ferrite.getfieldnames(sdh)
            ip = Ferrite.getfieldinterpolation(sdh, name)
            push!(fields, name=>ip)
        end
        push!(subdofhandlers, SerializedSubDofHandler(getcellset(sdh), fields))
    end
    return SerializedDofHandler(dh.grid, subdofhandlers)
end
