using Fleck
using Documenter

makedocs(;
    modules=[Fleck],
    authors="Andrew Dolgert <adolgert@uw.edu>",
    repo="https://github.com/adolgert/Fleck.jl/blob/{commit}{path}#L{line}",
    sitename="Fleck.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://adolgert.github.io/Fleck.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Structure" => "objects.md",
        "Develop" => "develop.md",
        "Delay Equations" => "delay.md",
        "Vector Addition Systems" => "vas.md",
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/adolgert/Fleck.jl",
)
