using NaturalNeighbours
using Documenter

DocMeta.setdocmeta!(NaturalNeighbours, :DocTestSetup, :(using NaturalNeighbours); recursive=true)

makedocs(;
    modules=[NaturalNeighbours],
    authors="Daniel VandenHeuvel <danj.vandenheuvel@gmail.com>",
    repo="https://github.com/DanielVandH/NaturalNeighbours.jl/blob/{commit}{path}#{line}",
    sitename="NaturalNeighbours.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DanielVandH.github.io/NaturalNeighbours.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Mathematical Details" => [
            "Interpolation" => "interpolation_math.md",
            "Differentiation" => "differentiation_math.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/DanielVandH/NaturalNeighbours.jl",
    devbranch="main",
)
