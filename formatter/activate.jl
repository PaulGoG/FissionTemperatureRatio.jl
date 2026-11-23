using Pkg
Pkg.activate(@__DIR__; io = devnull)
# The formatter is pinned to one release, so that the local gate and the CI gate agree. A manifest
# resolved under an earlier pin lies outside a raised one: `instantiate` does not notice it and
# `resolve` refuses it, so `update` moves it to the pinned release.
try
    Pkg.resolve(; io = devnull)
catch exception
    exception isa Pkg.Resolve.ResolverError || rethrow()
    Pkg.update(; io = devnull)
end
Pkg.instantiate(; io = devnull)
