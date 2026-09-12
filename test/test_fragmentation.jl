
@testset "fragmentation domain" begin
    A₀, Z₀ = 252, 98
    domain = fragmentation_domain(A₀, Z₀, 126:140, 5, FLAT_CHARGES)

    @testset "charge distribution is normalized per mass number" begin
        for A in unique(domain.A)
            @test sum(domain.p[domain.A .== A]) ≈ 1 atol = 1e-12
        end
    end

    @testset "both fragments of every pair are present" begin
        for A_H in 127:140
            @test A_H in domain.A
            @test A₀ - A_H in domain.A
        end
        @test length(charges(domain, 130)) == 5
    end

    @testset "the symmetric split is entered once" begin
        # Its complement is itself; entering it twice would double its weight.
        @test count(==(126), domain.A) == 5
    end

    @testset "charge polarization vanishes at the symmetric split" begin
        # The two fragments are the same nuclide there, so there is nothing to polarize, and the
        # retained charges must be invariant under Z -> Z₀ - Z.
        @test symmetric_charge_set_is_invariant(domain, A₀, Z₀) === true
        centre = most_probable_charge(126, A₀, Z₀, 0.0)
        @test centre ≈ Z₀ / 2
    end

    @testset "odd fissioning nucleus has no symmetric split" begin
        odd_domain = fragmentation_domain(235, 92, 118:130, 5, FLAT_CHARGES)
        @test symmetric_charge_set_is_invariant(odd_domain, 235, 92) === missing
    end

    @testset "invalid arguments are rejected" begin
        @test_throws ArgumentError fragmentation_domain(A₀, Z₀, 126:140, 4, FLAT_CHARGES)
        @test_throws ArgumentError fragmentation_domain(A₀, Z₀, 120:140, 5, FLAT_CHARGES)
        @test_throws ArgumentError fragmentation_domain(A₀, Z₀, 126:252, 5, FLAT_CHARGES)
    end

    @testset "averaging over charge skips missing values" begin
        values = Dict{Int,Union{Float64,Missing}}(
            Z => (Z == 49 ? missing : 2.0) for Z in charges(domain, 130)
        )
        @test average_over_charge(values, domain, 130) ≈ 2.0
        absent = Dict{Int,Union{Float64,Missing}}(Z => missing for Z in charges(domain, 130))
        @test ismissing(average_over_charge(absent, domain, 130))
    end
end
