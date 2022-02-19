using Fleck
using Documenter
using Literate

# Set Literate.jl config if not being compiled on recognized service.
config = Dict{String,String}()
if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
  config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/AlgebraicPetri.jl/blob/gh-pages/dev"
  config["repo_root_url"] = "https://github.com/AlgebraicJulia/AlgebraicPetri.jl/blob/master/docs"
end

const literate_dir = joinpath(@__DIR__, "..", "examples")
const generated_dir = joinpath(@__DIR__, "src", "examples")

for (root, dirs, files) in walkdir(literate_dir)
  out_dir = joinpath(generated_dir, relpath(root, literate_dir))
  for file in files
    f,l = splitext(file)
    if l == ".jl" && !startswith(f, "_")
      Literate.markdown(joinpath(root, file), out_dir;
        config=config, documenter=true, credit=false)
    end
  end
end

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
        "Develop" => "develop.md"
    ],
)

deploydocs(;
    repo="github.com/adolgert/Fleck.jl",
)
