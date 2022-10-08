using FerriteProblems
using Documenter

DocMeta.setdocmeta!(FerriteProblems, :DocTestSetup, :(using FerriteProblems); recursive=true)

makedocs(;
    modules=[FerriteProblems],
    authors="Knut Andreas Meyer and contributors",
    repo="https://github.com/KnutAM/FerriteProblems.jl/blob/{commit}{path}#{line}",
    sitename="FerriteProblems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://KnutAM.github.io/FerriteProblems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/KnutAM/FerriteProblems.jl",
    devbranch="main",
)
