include(joinpath(@__DIR__, "activate.jl"))

using Documenter
using DocumenterCitations
using FissionTemperatureRatio

DocMeta.setdocmeta!(
    FissionTemperatureRatio, :DocTestSetup, :(using FissionTemperatureRatio); recursive = true
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :numeric)

makedocs(;
    modules = [FissionTemperatureRatio],
    sitename = "FissionTemperatureRatio.jl",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    plugins = [bib],
    pages = [
        "Home" => "index.md",
        "Method" => "method.md",
        "Naming" => "naming.md",
        "Reference" => "reference.md",
        "References" => "references.md",
    ],
    checkdocs = :exports,
    # Stated rather than inferred from the git remote, so that the build also works from a source
    # tree without one.
    repo = Documenter.Remotes.GitHub("PaulGoG", "FissionTemperatureRatio.jl"),
)

deploydocs(; repo = "github.com/PaulGoG/FissionTemperatureRatio.jl", devbranch = "main")
