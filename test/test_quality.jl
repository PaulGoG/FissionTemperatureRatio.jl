using Aqua
using ExplicitImports
using JET
using CairoMakie: CairoMakie

@testset "package quality" begin
    Aqua.test_all(FissionTemperatureRatio; ambiguities = false)
    # The analysis covers the extension, which takes CairoMakie wholesale; every other dependency
    # is imported by name.
    @test isnothing(
        check_no_implicit_imports(FissionTemperatureRatio; skip = (Base, Core, CairoMakie))
    )
    @test isnothing(check_no_stale_explicit_imports(FissionTemperatureRatio))
    @test isnothing(check_all_qualified_accesses_via_owners(FissionTemperatureRatio))
    # Restricted to frames defined in this package: an unrestricted analysis walks into the
    # internals of Makie, CSV and the standard library, where the reports are not actionable here.
    JET.test_package(FissionTemperatureRatio; target_modules = (FissionTemperatureRatio,))

    extension = Base.get_extension(
        FissionTemperatureRatio, :FissionTemperatureRatioCairoMakieExt
    )
    @test extension !== nothing
    @test isnothing(check_no_stale_explicit_imports(extension))

    @testset "the FissionFragmentsDomain pin agrees across environments" begin
        root = pkgdir(FissionTemperatureRatio)
        project = TOML.parsefile(joinpath(root, "Project.toml"))
        rev = project["sources"]["FissionFragmentsDomain"]["rev"]
        compat = project["compat"]["FissionFragmentsDomain"]
        for d in ("test", "scripts", "docs", "bench")
            auxiliary = TOML.parsefile(joinpath(root, d, "Project.toml"))
            @test auxiliary["sources"]["FissionFragmentsDomain"]["rev"] == rev
            @test auxiliary["compat"]["FissionFragmentsDomain"] == compat
        end
    end
end
