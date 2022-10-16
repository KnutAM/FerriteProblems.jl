struct FerriteIO
    folder::ScalarWrapper{String}       # Where the results are saved
    def_file::String                    # Name of file saving the `def`
    postfile::String                    # Name of file saving the `post` (if given)
    nsteps_per_file::Int                # Maximum number of steps per data file
    switchsize::Float64                 # Size of data file (MB) when go to next file
    currentsize::ScalarWrapper{Int}     # Approximate current size of data file (B)
    time::Vector                        # Time history
    filesteps::Vector{Int}              # Step of datafiles step[k] gives first step of datafiles[k]
    datafiles::Vector{String}           # Name of each datafile
    fileobject::ScalarWrapper           # Current file, empty upon saving
    filenumber::ScalarWrapper{Int}      # Current file in `fileobject`
    dof2node::Dict{Symbol, Matrix{Int}} # Output from `reshape_to_nodes(1:ndof)`
end

"""
    FerriteIO(
        folder::String, def::FEDefinition, post=nothing; 
        def_file="FEDefinition.jld2", 
        postfile="FEPost.jld2",
        T=Float64, 
        nsteps_per_file=typemax(Int), 
        switchsize=Inf
        )

Constructor for creating a `FerriteIO` when simulating. 
"""
function FerriteIO(
    folder::String, def::FEDefinition, post=nothing; 
    def_file="FEDefinition.jld2", 
    postfile="FEPost.jld2",
    T=Float64, 
    nsteps_per_file=typemax(Int), 
    switchsize=Inf,
    )
    time = T[]
    firststep = [1,]
    datafiles=String[]
    dof2node = get_dof2node(def)
    isdir(folder) || mkdir(folder)
    jldsave(joinpath(folder, def_file); def=def)
    _postfile = isnothing(post) ? "" : postfile
    
    file = new_file!(firststep, datafiles, folder)
    return FerriteIO(
        ScalarWrapper(folder), def_file, _postfile, 
        nsteps_per_file, switchsize, ScalarWrapper(0),
        time, firststep, datafiles, 
        ScalarWrapper(file), ScalarWrapper(1),
        dof2node
        )
end

"""
    FerriteIO(filename::String)

Constructor for reading a `FerriteIO` that was saved during a simulation
"""
function FerriteIO(filename::String)
    key="io"
    endswith(filename, ".jld2") || @warn("Expected file ending with \".jld2\", unlike $filename")
    contents = load(filename)
    local io
    try
        io = contents[key]
    catch e
        isa(e, KeyError) ? KeyError("Expected $key to exist in $filename") : rethrow(e)
    end
    if !isa(io, FerriteIO)
        throw(ArgumentError("Expected $filename to have a FerriteIO object in the key \"$key\""))
    end
    io.folder[] = dirname(filename) # Update folder location, assume that other files maintain relative paths
    return io
end

"""
    filepath(io::FerriteIO, args...) = joinpath(io.folder[], args...)

Get the path of a file relative `io`'s folder
"""
filepath(io::FerriteIO, args...) = joinpath(io.folder[], args...)

"""
    datafilepath(io::FerriteIO, num=length(io.datafiles))

Get the path of the data file number `num` in `io.datafiles`
"""
datafilepath(io::FerriteIO, num=length(io.datafiles)) = filepath(io, io.datafiles[num])

"""
    close_problem(io::FerriteIO, post=nothing)

Method for closing all open files before ending the simulation
"""
function close_problem(io::FerriteIO, post=nothing)
    close(io.fileobject[])
    if !isnothing(post)
        jldsave(filepath(io, io.postfile); post=post)
    end
    jldsave(joinpath(io.folder[], "FerriteIO.jld2"), io=io)
end

function new_file!(io::FerriteIO)
    close(io.fileobject[])
    io.currentsize[] = 0
    io.fileobject[] = new_file!(io.filesteps, io.datafiles, io.folder[])
    io.filenumber[] += 1
end

function new_file!(filesteps::Vector{Int}, datafiles::Vector{String}, folder::String)
    push!(filesteps, filesteps[end])
    name = @sprintf("datafile_%05u", length(datafiles)+1)
    push!(datafiles, name)
    return jldopen(joinpath(folder, name), "w")
end

function new_file_if_needed!(io::FerriteIO)
    num_current_steps = io.filesteps[end]-io.filesteps[1]+1
    if num_current_steps > io.nsteps_per_file
        new_file!(io::FerriteIO)
    elseif filesize(datafilepath(io)) > io.switchsize*10^6  # switchsize given in MB 
        new_file!(io::FerriteIO)
    end
    return nothing
end

function addstep!(io::FerriteIO, time)
    io.filesteps[end] += 1   # Note: Must be done before new_file!
    new_file_if_needed!(io)  # Check size and number of steps
    push!(io.time, time)
    return nothing
end

function update_currentsize!(io, val)
    io.currentsize[] += Base.summarysize(val)
end

function savedata!(io::FerriteIO, vals, type::String, field::String, dt_order::Int)
    step = string(length(io.time))
    update_currentsize!(io, vals)
    file = io.fileobject[]
    file["$step/$dt_order/$type/$field"] = vals
end

const _dofkey = "dof"
const _nodekey = "node"
const _cellkey = "cell"
const _ipkey = "ip"
const _globalkey = "global"

function savedofdata!(io::FerriteIO, vals, dt_order=0)
    return savedata!(io, vals, _dofkey, _dofkey, dt_order)
end
function savenodedata!(io::FerriteIO, vals, field, dt_order=0)
    return savedata!(io, vals, _nodekey, field, dt_order)
end
function savecelldata!(io::FerriteIO, vals, field, dt_order=0)
    return savedata!(io, vals, _cellkey, field, dt_order)
end
function saveipdata!(io::FerriteIO, vals, field, dt_order=0)
    return savedata!(io, vals, _ipkey, field, dt_order)
end
function saveglobaldata!(io::FerriteIO, vals, field, dt_order=0)
    return savedata!(io, vals, _globalkey, field, dt_order)
end

getfilenumber(io, step) = findfirst(x -> x >= step, io.firststep)

function open_if_needed!(io::FerriteIO, step, mode="r")
    num = getfilenumber(io, step)
    if num != io.filenumber[]
        close(io.fileobject[])
        io.filenumber[] = jldopen(datafilepath(io, num), mode)
        io.filenumber[] = num
    end
end

function checkkey(file, step, dt_order, type, field)
    # Do full check first, to avoid overhead of going through each level
    haskey(file, "$step/$dt_order/$type/$field") && return nothing
    if !haskey(file, "$step")
        throw(ArgumentError("Step $step is not saved"))
    end
    stepgroup = file["$step"]
    if !haskey(stepgroup, "$dt_order")
        throw(ArgumentError("Step $step has no $(dt_order)th order derivative data saved"))
    end
    ordergroup = stepgroup["$dt_order"]
    if !haskey(ordergroup, "$type")
        throw(ArgumentError("Step $step has no $(dt_type). order derivative $type-data saved"))
    end
    typegroup = ordergroup["$type"]
    if !haskey(typegroup, "$field")
        throw(ArgumentError("Step $step has not saved the $(dt_type). order derivative of $field of $type-data"))
    end
    throw(ErrorException("Unexpected to reach here"))
end

function getdata(io::FerriteIO, step::Int, type, field, dt_order)
    open_if_needed!(io::FerriteIO, step)
    file = io.fileobject[]
    checkkey(file, step, dt_order, type, field) # Errors if not ok
    return file["$step/$dt_order/$type/$field"]
end

function getdofdata(io::FerriteIO, step; dt_order=0)
    getdata(io, step, _dofkey, _dofkey, dt_order)
end
function getnodedata(io::FerriteIO, step, field; dt_order=0)
    getdata(io, step, _nodekey, field, dt_order)
end
function getcelldata(io::FerriteIO, step, field; dt_order=0)
    getdata(io, step, _cellkey, field, dt_order)
end
function getipdata(io::FerriteIO, step, field; dt_order=0)
    getdata(io, step, _ipkey, field, dt_order)
end
function getglobaldata(io::FerriteIO, step, field; dt_order=0)
    getdata(io, step, _globalkey, field, dt_order)
end


# Utility functions
get_dof2node(def::FEDefinition) = get_dof2node(getdh(def))
function get_dof2node(dh)
    fieldnames = Ferrite.getfieldnames(dh)
    dofs = collect(1:ndofs(dh))
    return Dict(field=>reshape_to_nodes(dh, dofs, field) for field in fieldnames)
end

