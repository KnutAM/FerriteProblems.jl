using FerriteProblems
using Documenter

const is_ci = get(ENV, "CI", "false") == "true"

include("generate.jl")
examples = ["plasticity.jl",]
GENERATEDEXAMPLES = [joinpath("examples", replace(f, ".jl"=>".md")) for f in examples]

build_examples(examples)

DocMeta.setdocmeta!(FerriteProblems, :DocTestSetup, :(using FerriteProblems); recursive=true)

makedocs(;
    modules=[FerriteProblems],
    authors="Knut Andreas Meyer and contributors",
    repo="https://github.com/KnutAM/FerriteProblems.jl/blob/{commit}{path}#{line}",
    sitename="FerriteProblems.jl",
    format=Documenter.HTML(;
        prettyurls=is_ci,
        canonical="https://KnutAM.github.io/FerriteProblems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => GENERATEDEXAMPLES,  
    ],
)

deploydocs(;
    repo="github.com/KnutAM/FerriteProblems.jl",
    devbranch="main",
    push_preview=true,
)
