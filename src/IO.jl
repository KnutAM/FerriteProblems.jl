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
    filesteps = [1,]
    datafiles=String[]
    dof2node = get_dof2node(def)
    isdir(folder) || mkdir(folder)
    jldsave(joinpath(folder, def_file); def=def)
    _postfile = isnothing(post) ? "" : postfile
    
    file = new_file!(filesteps, datafiles, folder)
    return FerriteIO(
        ScalarWrapper(folder), def_file, _postfile, 
        nsteps_per_file, switchsize, ScalarWrapper(0),
        time, filesteps, datafiles, 
        ScalarWrapper(file), ScalarWrapper(1),
        dof2node
        )
end

"""
    FerriteIO(filename::String)

Constructor for reading a `FerriteIO` that was saved during a simulation
"""
function FerriteIO(filename::String)
    endswith(filename, ".jld2") || @warn("Expected file ending with \".jld2\", unlike $filename")
    io = get_problem_part(filename, "io", FerriteIO)
    io.folder[] = dirname(filename) # Update folder location, assume that other files maintain relative paths
    io.filenumber[] = -1            # Force re-opening file if needed
    return io
end

"""
    filepath(io::FerriteIO, args...) = joinpath(io.folder[], args...)

Get the path of a file relative `io`'s folder
"""
filepath(io::FerriteIO, args...) = joinpath(io.folder[], args...) # internal

"""
    datafilepath(io::FerriteIO, num=length(io.datafiles))

Get the path of the data file number `num` in `io.datafiles`
"""
datafilepath(io::FerriteIO, num=length(io.datafiles)) = filepath(io, io.datafiles[num]) # internal

"""
    close_io(io::FerriteIO, post)

Close the currently open file in `io`, then the postprocessing 
struct to a jld2 file (if not `post != nothing`),
before finally saving the current `io` object to a .jld2 file
"""
function close_io(io::FerriteIO, post)   # internal
    close(io.fileobject[])
    isnothing(post) || jldsave(filepath(io, io.postfile); post=post)
    jldsave(joinpath(io.folder[], "FerriteIO.jld2"), io=io)
end
# Do nothing if there is no io defined. 
close_io(args...) = nothing 

"""
    new_file!
"""
function new_file!(io::FerriteIO)   # internal
    close(io.fileobject[])
    io.currentsize[] = 0
    io.fileobject[] = new_file!(io.filesteps, io.datafiles, io.folder[])
    io.filenumber[] += 1
end

function new_file!(filesteps::Vector{Int}, datafiles::Vector{String}, folder::String) # internal
    push!(filesteps, filesteps[end])
    name = @sprintf("datafile_%05u.jld2", length(datafiles)+1)
    push!(datafiles, name)
    return jldopen(joinpath(folder, name), "w")
end

"""
    new_file_if_needed!
"""
function new_file_if_needed!(io::FerriteIO) # internal
    num_current_steps = io.filesteps[end]-io.filesteps[1]+1
    if num_current_steps > io.nsteps_per_file
        new_file!(io::FerriteIO)
    elseif filesize(datafilepath(io)) > io.switchsize*10^6  # switchsize given in MB 
        new_file!(io::FerriteIO)
    end
    return nothing
end

"""
    FerriteProblems.addstep!(io::FerriteIO, p::FerriteProblem)

Add a new step to be saved by `io` at the time `gettime(p)`
Must be called before adding any new data
"""
function addstep!(io::FerriteIO, time)
    io.filesteps[end] += 1   # Note: Must be done before new_file!
    new_file_if_needed!(io)  # Check size and number of steps
    push!(io.time, time)
    return nothing
end

"""
    update_currentsize!
"""
function update_currentsize!(io, val) # internal
    io.currentsize[] += Base.summarysize(val)
end

"""
    savedata!
"""
function savedata!(io::FerriteIO, vals, type::String, field::String, dt_order::Int) # internal
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

"""
    FerriteProblems.savedofdata!(io::FerriteIO, vals, dt_order=0, field=\"$_dofkey\")

Save data pertaining to each degree of freedom.
Use a different field than `\"$_dofkey\"`` to save data located at each dof, 
but not the actual dof values (e.g. the residual vector)
"""
function savedofdata!(io::FerriteIO, vals, dt_order=0, field=_dofkey)
    return savedata!(io, vals, _dofkey, field, dt_order)
end

"""
    FerriteProblems.savenodedata!(io::FerriteIO, vals, field, dt_order=0)

Save data located at each node.
By convention this should be indexed by the node numbers in the grid.
(E.g. a `Vector` for all nodes or a `Dict{Int}` with keys the node numbers)
"""
function savenodedata!(io::FerriteIO, vals, field, dt_order=0)
    return savedata!(io, vals, _nodekey, field, dt_order)
end

"""
    FerriteProblems.savecelldata!(io::FerriteIO, vals, field, dt_order=0)

Save data for each cell. 
By convention this should be indexed by the cell numbers in the grid.
(E.g. a `Vector` for all cells or a `Dict{Int}` with keys the cell indices)
"""
function savecelldata!(io::FerriteIO, vals, field, dt_order=0)
    return savedata!(io, vals, _cellkey, field, dt_order)
end

"""
    FerriteProblems.saveipdata!(io::FerriteIO, vals, field, dt_order=0)

Save data for each integration point in cells in the grid. 
By convention the data for each cell should be indexed by the cell numbers in the grid.
(E.g. a `Vector` for all cells or a `Dict{Int}` with keys the cell indices)
Note that it is on the user to know how the integration points are numbered, 
i.e. which `QuadratureRule` that was used. 
"""
function saveipdata!(io::FerriteIO, vals, field, dt_order=0)
    return savedata!(io, vals, _ipkey, field, dt_order)
end

"""
    FerriteProblems.saveglobaldata!(io::FerriteIO, vals, field, dt_order=0)

Save data that is global to the entire simulation, i.e. global quantites such as 
reaction forces, total dissipation, etc. 
"""
function saveglobaldata!(io::FerriteIO, vals, field, dt_order=0)
    return savedata!(io, vals, _globalkey, field, dt_order)
end

"""
    getfilenumber
"""
getfilenumber(io, step) = findfirst(x -> x > step, io.filesteps) - 1    # internal

"""
    open_if_needed!
"""
function open_if_needed!(io::FerriteIO, step, mode="r") # internal
    num = getfilenumber(io, step)
    if num != io.filenumber[]
        close(io.fileobject[])
        io.fileobject[] = jldopen(datafilepath(io, num), mode)
        io.filenumber[] = num
    end
end

"""
    checkkey
"""
function checkkey(file, step, dt_order, type, field) # internal
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

"""
    gettimedata(io::FerriteIO)
    gettimedata(io::FerriteIO, step)
"""
gettimedata(io::FerriteIO) = io.time
gettimedata(io::FerriteIO, step) = io.time[step]

"""
    getdata
"""
function getdata(io::FerriteIO, step::Int, type, field, dt_order)
    open_if_needed!(io::FerriteIO, step)
    file = io.fileobject[]
    checkkey(file, step, dt_order, type, field) # Errors if not ok
    # Copy, otherwise we can change the contents!!!
    return copy(file["$step/$dt_order/$type/$field"])
end

"""
    FerriteProblems.getdofdata(io::FerriteIO, step, field=\"$_dofkey\"; dt_order=0)

Get the data saved by [`FerriteProblems.savedofdata!`](@ref) in `step` for `field`.
"""
function getdofdata(io::FerriteIO, step; dt_order=0, field=_dofkey)
    getdata(io, step, _dofkey, field, dt_order)
end

"""
    FerriteProblems.getnodedata(io::FerriteIO, step, field; dt_order=0)

Get the data saved by [`FerriteProblems.savenodedata!`](@ref) in `step` for `field`.
"""
function getnodedata(io::FerriteIO, step, field; dt_order=0)
    getdata(io, step, _nodekey, field, dt_order)
end

"""
    FerriteProblems.getcelldata(io::FerriteIO, step, field; dt_order=0)

Get the data saved by [`FerriteProblems.savecelldata!`](@ref) in `step` for `field`.
"""
function getcelldata(io::FerriteIO, step, field; dt_order=0)
    getdata(io, step, _cellkey, field, dt_order)
end

"""
    FerriteProblems.getipdata(io::FerriteIO, step, field; dt_order=0)

Get the data saved by [`FerriteProblems.saveipdata!`](@ref) in `step` for `field`.
"""
function getipdata(io::FerriteIO, step, field; dt_order=0)
    getdata(io, step, _ipkey, field, dt_order)
end

"""
    FerriteProblems.getglobaldata(io::FerriteIO, step, field; dt_order=0)

Get the data saved by [`FerriteProblems.saveglobaldata!`](@ref) in `step` for `field`.
"""
function getglobaldata(io::FerriteIO, step, field; dt_order=0)
    getdata(io, step, _globalkey, field, dt_order)
end

function get_problem_part(filename, key, Type)
    contents = load(filename)
    local part
    try
        part = contents[key]
    catch e
        isa(e, KeyError) ? KeyError("Expected $key to exist in $filename") : rethrow(e)
    end
    if !isa(part, Type)
        throw(ArgumentError("Expected $filename to have a $Type object in the key \"$key\""))
    end
    return part
end

"""
    FerriteProblems.getdef(io::FerriteIO)

Load the `FEDefinition` from the results saved by `io`
"""
function getdef(io::FerriteIO)
    filename = joinpath(io.folder[], io.def_file)
    return get_problem_part(filename, "def", FEDefinition)
end

"""
    FerriteProblems.getpost(io::FerriteIO)

Load the user defined `post` from the results saved by `io`
"""
function getpost(io::FerriteIO)
    filename = joinpath(io.folder[], io.postfile)
    return get_problem_part(filename, "post", Any)
end
    

# Utility functions
"""
    get_dof2node
"""
get_dof2node(def::FEDefinition) = get_dof2node(getdh(def))
function get_dof2node(dh)
    fieldnames = Ferrite.getfieldnames(dh)
    dofs = collect(1:ndofs(dh))
    return Dict(field=>reshape_to_nodes(dh, dofs, field) for field in fieldnames)
end

