# generate examples
import Literate

# Recursive reading for deleting output files
function readdir_recursive(args...)
    items = readdir(args...; join=true)
    subitems = eltype(items)[]
    foreach(item->isdir(item) ? append!(subitems, readdir_recursive(item)) : nothing, items)
    return append!(items, subitems)
end
function remove_generalted_results(args...)
    should_delete(file) = any(endswith(file, ext) for ext in args)
    for dirname in ("src", "build")
        thedir = joinpath(@__DIR__, dirname, "examples")
        foreach(file ->  should_delete(file) && rm(file), readdir_recursive(thedir))
    end
end

function build_examples(examples)
    EXAMPLEDIR = joinpath(@__DIR__, "src", "literate")
    GENERATEDDIR = joinpath(@__DIR__, "src", "examples")
    rm(GENERATEDDIR; force=true, recursive=true)
    mkpath(GENERATEDDIR)

    # Copy supplementary files first
    suplementary_fileextensions = [".inp", ".svg", ".png", ".jpg", ".gif"]
    for example in readdir(EXAMPLEDIR)
        supplementary_ext = any(endswith.(example, suplementary_fileextensions))
        supplementary_jl = endswith(example, ".jl") && example ∉ examples
        if supplementary_ext || supplementary_jl || !isfile(example)
            cp(joinpath(EXAMPLEDIR, example), joinpath(GENERATEDDIR, example); force=true)
        end
    end

    for example in examples
        input = abspath(joinpath(EXAMPLEDIR, example))
        isfile(input) || throw(SystemError("$input not found"))
        script = Literate.script(input, GENERATEDDIR)
        code = strip(read(script, String))

        # remove "hidden" lines which are not shown in the markdown
        line_ending_symbol = occursin(code, "\r\n") ? "\r\n" : "\n"
        code_clean = join(filter(x->!endswith(x,"#hide"),split(code, r"\n|\r\n")), line_ending_symbol)

        mdpost(str) = replace(str, "@__CODE__" => code_clean)
        Literate.markdown(input, GENERATEDDIR, postprocess = mdpost)
        Literate.notebook(input, GENERATEDDIR, execute = is_ci) # Don't execute locally
    end

end
