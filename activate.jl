using Pkg
# FissionFragmentsDomain enters by URL. The git executable honours the user's git configuration
# (credentials, `url.insteadOf` rewrites to SSH), where libgit2 can hang on the fetch.
haskey(ENV, "JULIA_PKG_USE_CLI_GIT") || (ENV["JULIA_PKG_USE_CLI_GIT"] = "true")
Pkg.activate(@__DIR__; io = devnull)
# A changed `[sources]` revision leaves an existing manifest stale; `instantiate` alone does not
# notice.
Pkg.resolve(; io = devnull)
Pkg.instantiate(; io = devnull)
