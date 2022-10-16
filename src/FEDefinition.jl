struct FEDefinition{DH,CH,NH,CV,M,ST}
    # FE-handlers
    dh::DH  # DofHandler/MixedDofHandler
    ch::CH  # ConstraintHandler
    nh::NH  # NeumannHandler
    # Items related to each cell type in the grid
    # If multiple cell types, these should be tuples
    cv::CV  # cellvalues
    m::M    # material
    initialstate::ST
end

getdh(def::FEDefinition) = def.dh
getch(def::FEDefinition) = def.ch
getnh(def::FEDefinition) = def.nh
getcv(def::FEDefinition) = def.cv
getmaterial(def::FEDefinition) = def.m