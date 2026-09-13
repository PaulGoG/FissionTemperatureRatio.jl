using Pkg
Pkg.activate(@__DIR__; io = devnull)
# `[sources]` is honoured from Pkg 1.11 onwards. On the declared 1.10 floor it is silently
# discarded, and the unregistered package would then be looked up in a registry and not found, so
# below 1.11 it is developed by path instead.
if VERSION < v"1.11"
    Pkg.develop(; path = dirname(@__DIR__), io = devnull)
end
Pkg.instantiate(; io = devnull)
