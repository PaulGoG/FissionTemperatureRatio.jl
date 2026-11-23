using Pkg
# FissionFragmentsDomain enters by URL. The git executable honours the user's git configuration
# (credentials, `url.insteadOf` rewrites to SSH), where libgit2 can hang on the fetch.
haskey(ENV, "JULIA_PKG_USE_CLI_GIT") || (ENV["JULIA_PKG_USE_CLI_GIT"] = "true")
Pkg.activate(@__DIR__; io = devnull)
# A new depot holds no registry, and `resolve` does not fetch one.
isempty(Pkg.Registry.reachable_registries()) && Pkg.Registry.add()
# The package enters this environment by path, so a dependency added to it leaves the resolved
# manifest here stale; `instantiate` alone does not notice.
Pkg.resolve(; io = devnull)
Pkg.instantiate(; io = devnull)
