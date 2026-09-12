include(joinpath(@__DIR__, "activate.jl"))

using Documenter
using FissionTemperatureRatio

DocMeta.setdocmeta!(
    FissionTemperatureRatio, :DocTestSetup, :(using FissionTemperatureRatio); recursive = true
)

makedocs(;
    modules = [FissionTemperatureRatio],
    sitename = "FissionTemperatureRatio.jl",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    pages = ["Home" => "index.md", "Method" => "method.md", "Reference" => "reference.md"],
    checkdocs = :exports,
    # The repository is local, so there is no remote to link source lines against.
    remotes = nothing,
)
