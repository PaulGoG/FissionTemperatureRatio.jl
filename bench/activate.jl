using Pkg
Pkg.activate(@__DIR__; io = devnull)
# The package enters this environment by path, so a dependency added to it leaves the resolved
# manifest here stale; `instantiate` alone does not notice.
Pkg.resolve(; io = devnull)
Pkg.instantiate(; io = devnull)
