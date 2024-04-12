using Fleck
using Documenter
using Literate


rebuild_literate = "--fresh" in ARGS

# Literate puts images into png files in the same directory as the source
# but we would like them in the assets subdirectory instead, so here we
# postprocess the markdown generated by Literate to tell it where to look.
function postprocess_markdown(sfile::String)
    replace(sfile, r"\(([a-z0-9_\-]+\.(png|pdf|svg))\)" => s"(literate/\1)")
end


example_base = joinpath(dirname(@__FILE__), "src")
adliterate = [
        ("simple_board.jl", "mainloop"),
        ("distributions.jl", "distributions"), 
        ("constant_birth.jl", "constant_birth"),
        ("sir.jl", "sir"),
        ("commonrandom.jl", "commonrandom")
    ]
literate_subdir = joinpath(example_base, "literate")
isdir(literate_subdir) || mkdir(literate_subdir)

if rebuild_literate
    for rmchunk in readdir(literate_subdir; join=true)
        rm(rmchunk)
    end
end

for (source, target) in adliterate
    fsource, ftarget = joinpath.(example_base, [source, target])
    if rebuild_literate
        rm("$(ftarget).md", force=true)
    end
    if !isfile("$(ftarget).md") || mtime(fsource) > mtime("$(ftarget).md")
        println("Literate of $source to $target")
        Literate.markdown(
            fsource,
            example_base,
            name=target,
            postprocess=postprocess_markdown,
            execute=true
            )
        ischunk = Regex("$(target)-[0-9]+.(png|svg|pdf)")
        chunks = [fn for fn in readdir(example_base) if match(ischunk, fn) !== nothing]
        for chunk in chunks
            mv(joinpath(example_base, chunk), joinpath(example_base, "literate", chunk))
        end
    end
end

makedocs(;
    modules=[Fleck],
    authors="Andrew Dolgert <adolgert@uw.edu>",
    sitename="Fleck.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        size_threshold_warn=2^17,
        canonical="https://adolgert.github.io/Fleck.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => [
            "guide.md",
            "mainloop.md",
            "distributions.md"
        ],
        "Manual" => [
            "Structure" => "objects.md",
            "commonrandom.md",
            "background.md",
            "distrib.md",
            "Develop" => "develop.md",
            "samplers.md",
            "Vector Addition Systems" => "vas.md",
        ],
        "Examples" => [
            "Birth-death Process" => "constant_birth.md",
            "SIR Model" => "sir.md"
        ],
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    target = "build",
    repo = "github.com/adolgert/Fleck.jl.git",
    branch = "gh-pages"
)
