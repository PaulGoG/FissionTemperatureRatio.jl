
@testset "fragmentation domain" begin
    A₀, Z₀ = 252, 98
    domain = fragmentation_domain(A₀, Z₀, 126:140, 5, FLAT_CHARGES)

    @testset "charge distribution is the analytic Gaussian, not renormalized" begin
        # Each weight is the Gaussian as evaluated, carrying its own normalizing factor. This is
        # the invariant: the distribution is not rescaled to sum to one over the charges retained.
        Zₚ = most_probable_charge(130, A₀, Z₀, FLAT_CHARGES.fallback_ΔZ)
        rms = FLAT_CHARGES.fallback_rms
        for Z in charges(domain, 130)
            @test charge_probability(domain, 130, Z) ≈
                exp(-(Z - Zₚ)^2 / (2 * rms^2)) / (sqrt(2π) * rms)
        end

        # Summed over a mass number the result is near unity but not equal to it, two effects
        # pulling opposite ways: charges outside the retained window are lost, while evaluating a
        # density on a unit charge lattice counts the peak more heavily than integrating it. At
        # five charges and rms 0.6 the window spans more than three dispersions, so the lattice
        # term dominates and the sum sits just above one.
        for A in unique(domain.A)
            @test sum(domain.p[domain.A .== A]) ≈ 1 atol = 5.0e-3
            @test sum(domain.p[domain.A .== A]) != 1
        end

        # Retaining too few charges per mass number is meant to cost something visible, which is
        # the reason for not renormalizing: one charge per mass spans a third of the distribution.
        narrow = fragmentation_domain(A₀, Z₀, 126:140, 1, FLAT_CHARGES)
        @test sum(narrow.p[narrow.A .== 130]) < 0.7
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
