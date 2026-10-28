using Pkg
Pkg.activate(@__DIR__; io = devnull)
# The formatter is pinned to one release, so that the local gate and the CI gate agree; a changed
# pin leaves an existing manifest stale, which `instantiate` alone does not notice.
Pkg.resolve(; io = devnull)
Pkg.instantiate(; io = devnull)
