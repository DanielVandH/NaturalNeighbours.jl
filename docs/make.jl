using NaturalNeighbourInterp
using Documenter

DocMeta.setdocmeta!(NaturalNeighbourInterp, :DocTestSetup, :(using NaturalNeighbourInterp); recursive=true)

makedocs(;
    modules=[NaturalNeighbourInterp],
    authors="Daniel VandenHeuvel <danj.vandenheuvel@gmail.com>",
    repo="https://github.com/DanielVandH/NaturalNeighbourInterp.jl/blob/{commit}{path}#{line}",
    sitename="NaturalNeighbourInterp.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DanielVandH.github.io/NaturalNeighbourInterp.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DanielVandH/NaturalNeighbourInterp.jl",
    devbranch="main",
)
