using Aqua
using ExplicitImports
using JET
using CairoMakie: CairoMakie

@testset "package quality" begin
    Aqua.test_all(FissionTemperatureRatio; ambiguities = false)
    # Makie is designed to be brought in wholesale; every other dependency is imported by name.
    @test isnothing(
        check_no_implicit_imports(FissionTemperatureRatio; skip = (Base, Core, CairoMakie))
    )
    @test isnothing(check_no_stale_explicit_imports(FissionTemperatureRatio))
    @test isnothing(check_all_qualified_accesses_via_owners(FissionTemperatureRatio))
    # Restricted to frames defined in this package: an unrestricted analysis walks into the
    # internals of Makie, CSV and the standard library, where the reports are not actionable here.
    JET.test_package(FissionTemperatureRatio; target_modules = (FissionTemperatureRatio,))
end
