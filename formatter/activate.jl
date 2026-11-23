using Pkg
Pkg.activate(@__DIR__; io = devnull)
# A new depot holds no registry, and `resolve` does not fetch one.
isempty(Pkg.Registry.reachable_registries()) && Pkg.Registry.add()
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
