using FerriteProblems
using Documenter

const is_ci = get(ENV, "CI", "false") == "true"

include("generate.jl")
examples = ["plasticity.jl", "io.jl", "transient_heat.jl", "incompressible_elasticity.jl", "porous_media.jl"]
GENERATEDEXAMPLES = [joinpath("examples", replace(f, ".jl"=>".md")) for f in examples]

# Build examples, see `generate.jl`
build_examples(examples)

DocMeta.setdocmeta!(FerriteProblems, :DocTestSetup, :(using FerriteProblems); recursive=true)

#= Temp testing custom macros...
jax = Documenter.MathJax(Dict(:TeX => Dict(
    :equationNumbers => Dict(:autoNumber => "AMS"),
    :Macros => Dict(
        :intO => ["{\\int_\\Omega #1 \\ \\mathrm{d}\\Omega}", 1],
        :mytest => "a"
    ),
)))

katex=Documenter.KaTeX(Dict(:Macros=>Dict(
    :intO => ["{\\int_\\Omega #1 \\ \\mathrm{d}\\Omega}", 1],
    :mytest => "a"
    )))
=#

makedocs(;
    authors="Knut Andreas Meyer and contributors",
    sitename="FerriteProblems.jl",
    format=Documenter.HTML(;
        canonical="https://KnutAM.github.io/FerriteProblems.jl",
        assets=String[],
        collapselevel = 1,
    ),
    warnonly = true,
    doctest = false,
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "Examples" => GENERATEDEXAMPLES,
        "Internals" => "internals.md",
    ],
)

# Remove any generated files, see `generate.jl`
remove_generalted_results(".vtu", ".pvd", ".jld2")

deploydocs(;
    repo="github.com/KnutAM/FerriteProblems.jl",
    devbranch="main",
    push_preview=true,
)
