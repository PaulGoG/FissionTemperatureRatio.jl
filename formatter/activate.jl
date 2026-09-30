using Pkg
# Git fetches through the git executable, as in every other environment of the repository.
haskey(ENV, "JULIA_PKG_USE_CLI_GIT") || (ENV["JULIA_PKG_USE_CLI_GIT"] = "true")
Pkg.activate(@__DIR__; io = devnull)
Pkg.instantiate(; io = devnull)
