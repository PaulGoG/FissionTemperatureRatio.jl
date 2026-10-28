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

    @testset "write_figures writes one file per content key" begin
        # The synthetic run of the pipeline tests.
        mktempdir() do directory
            open(joinpath(directory, "shell_corrections.dat"), "w") do io
                println(io, "n S_N S_Z")
                for n in 20:119
                    println(io, n, " 0.0 0.0")
                end
                println(io, 120, " 1.0 1.0")
            end
            mkpath(joinpath(directory, "datasets"))
            write_multiplicity(joinpath(directory, "datasets", "full.dat"), 126:140)
            write_multiplicity(joinpath(directory, "datasets", "partial.dat"), 131:140)
            write_multiplicity(joinpath(directory, "datasets", "sparse.dat"), 126:4:138)
            path = joinpath(directory, "configuration.toml")
            write(path, PIPELINE_CONFIGURATION)
            result = run_pipeline(
                load_configuration(path; data_directory = directory); leave_one_out = false
            )

            figures = joinpath(directory, "figures")
            written = write_figures(result, figures)
            @test !isempty(written)
            @test all(isfile, values(written))
            @test all(endswith(".pdf"), values(written))
        end
    end
end
