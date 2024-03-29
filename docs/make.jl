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
        "Examples" => [
            "Interpolation" => "interpolation.md",
            "Differentiation" => "differentiation.md",
            "Switzerland Elevation Data" => "swiss.md"
        ],
        "Comparison of Interpolation Methods" => "compare.md",
        "Mathematical Details" => [
            "Interpolation Details" => "interpolation_math.md",
            "Differentiation Details" => "differentiation_math.md"
        ]
    ],
    warnonly=true
)

deploydocs(;
    repo="github.com/DanielVandH/NaturalNeighbours.jl",
    devbranch="main",
    push_preview=true
)
