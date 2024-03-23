using Fleck
using Documenter
using Literate

example_base = joinpath(dirname(@__FILE__), "src")
println("example base is $example_base")
Literate.markdown(
    joinpath(example_base, "simple_board.jl"),
    example_base,
    name="mainloop",
    execute=true
    )

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
        "Guide" => [
            "guide.md",
            "mainloop.md",
            "distributions.md",
            "rules.md"
        ],
        "Manual" => [
            "Structure" => "objects.md",
            "Delay Equations" => "delay.md",
            "Develop" => "develop.md",
            "Vector Addition Systems" => "vas.md",
        ],
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/adolgert/Fleck.jl",
)
