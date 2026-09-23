
@testset "multiplicity ratio" begin
    A₀ = 252
    data = Multiplicity([120, 126, 132], [1.0, 2.0, 3.0], [0.1, 0.2, 0.3], "synthetic", "-")

    @testset "definition and symmetry" begin
        ratio = multiplicity_ratio(data, A₀, 126:132)
        # The symmetric split pairs a fragment with itself, so the ratio is one half exactly.
        index = findfirst(==(126), ratio.A_H)
        @test ratio.ratio[index] ≈ 0.5
        # ν(132) = 3 against its complement ν(120) = 1.
        index = findfirst(==(132), ratio.A_H)
        @test ratio.ratio[index] ≈ 3 / 4
    end

    @testset "uncertainty propagation" begin
        ratio = multiplicity_ratio(data, A₀, 132:132)
        ν_H, σ_H, ν_L, σ_L = 3.0, 0.3, 1.0, 0.1
        expected = sqrt((ν_L * σ_H)^2 + (ν_H * σ_L)^2) / (ν_L + ν_H)^2
        @test only(ratio.σ) ≈ expected
    end

    @testset "unpaired and unphysical points are omitted, never clipped" begin
        unpaired = Multiplicity([132], [3.0], [missing], "unpaired", "-")
        @test isempty(multiplicity_ratio(unpaired, A₀, 126:140))
        # A vanishing light multiplicity puts the ratio at the singular endpoint of the
        # temperature ratio relation; the point must be dropped rather than pushed to 1.
        degenerate = Multiplicity([120, 132], [0.0, 3.0], [missing, missing], "degenerate", "-")
        ratio = multiplicity_ratio(degenerate, A₀, 132:132)
        @test isempty(ratio)
    end

    @testset "reader stores an absent uncertainty as missing, keeping the point" begin
        mktempdir() do directory
            path = joinpath(directory, "set.dat")
            write(path, "A nu nu_uncertainty\n126 2.0 0.1\n132 3.0\n")
            data = read_multiplicity(path; label = "set")
            @test length(data) == 2
            @test ismissing(data.σν[2])
        end
    end
end
