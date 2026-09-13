@testset "plotting extension" begin
    extension_name = :FissionTemperatureRatioCairoMakieExt

    # The plotting stack must stay weak: this is what keeps `using FissionTemperatureRatio` free
    # of graphics. Asserted on the manifest of the package rather than on load order, because the
    # quality checks earlier in this suite already load extensions of their own accord.
    project = TOML.parsefile(joinpath(pkgdir(FissionTemperatureRatio), "Project.toml"))
    @test !haskey(project["deps"], "CairoMakie")
    @test haskey(project["weakdeps"], "CairoMakie")
    @test haskey(project["extensions"], String(extension_name))

    @eval using CairoMakie
    extension = Base.get_extension(FissionTemperatureRatio, extension_name)
    @test extension !== nothing

    for stub in (publication_theme, plot_multiplicities, plot_ratio, save_figure, write_figures)
        @test !isempty(methods(stub))
    end

    @testset "figures are produced and written" begin
        data = Multiplicity([120, 132], [3.10, 0.69], [0.05, 0.03], "example", "")
        curve = RatioCurve([130, 132], [1.20, 1.10], [0.02, 0.01], "example")

        CairoMakie.with_theme(publication_theme()) do
            figure = plot_multiplicities([data]; A_0 = 252)
            @test figure isa CairoMakie.Figure

            ratio = plot_ratio(
                [curve]; ylabel = "R", reference = 1.0, reference_label = "Unity"
            )
            @test ratio isa CairoMakie.Figure

            mktempdir() do directory
                path = joinpath(directory, "nested", "figure.pdf")
                written = save_figure(path, ratio)
                @test isfile(written)
                @test filesize(written) > 0
                # An existing file is never overwritten; a suffixed path is used instead.
                again = save_figure(path, ratio)
                @test again != written
                @test isfile(again)
            end
            return nothing
        end
    end
end
