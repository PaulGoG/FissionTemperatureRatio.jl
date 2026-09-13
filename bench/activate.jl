using Pkg
Pkg.activate(@__DIR__; io = devnull)
# Developed by path rather than declared in `[sources]`: that key is honoured only from Pkg 1.11,
# and on the declared floor `Pkg.test` refuses to merge a test project that carries it.
Pkg.develop(; path = dirname(@__DIR__), io = devnull)
Pkg.instantiate(; io = devnull)
