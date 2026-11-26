# The extraction inverts FissionFragmentsDomain's relation between R_T and E*_H/TXE, which the
# consuming emission code applies forwards. These tests hold the two directions together on the
# shipped domains; the relation itself is tested in that package.

"The multiplicity ratio a partition with temperature ratio `R_T` gives at every heavy mass."
function forward_ratio(averaging, model, domain, R_T)
    A_H = Int[]
    r = Float64[]
    for mass in domain.heavy_masses
        value = heavy_excitation_fraction(averaging, model, domain, mass, R_T)
        value === nothing && continue
        push!(A_H, mass)
        push!(r, value)
    end
    return RatioCurve(A_H, r, fill(missing, length(A_H)), "forward")
end

@testset "temperature ratio" begin
    @testset "the extraction inverts the shared forward relation to rounding" begin
        # For every heavy mass of the 252-Cf and 235-U domains the extraction, applied to the
        # E*_H/TXE a partition with temperature ratio R_T gives, returns R_T: the round trip a
        # consuming code closes when it partitions every fragment pair with its own a_L/a_H.
        for label in ("Cf252_sf", "U235_nth")
            domain = shipped_domain(label)
            masses = TEST_MASSES
            weighted = ChargeResolved(
                mean_total_excitation(masses, domain, synthetic_kinetic_energy(domain))
            )
            for model in (BSFG_MODEL, GC_MODEL),
                averaging in (ChargeResolved(), weighted),
                R_T in (0.8, 1.0, 1.45)

                forward = forward_ratio(averaging, model, domain, R_T)
                @test forward.A_H == collect(domain.heavy_masses)
                back = temperature_ratio(averaging, model, domain, forward)
                @test back.A_H == forward.A_H
                @test all(isapprox.(back.ratio, R_T; rtol = 1e-12))
            end
        end
    end

    @testset "r_ν = 1/2 at the symmetric split gives R_T = 1" begin
        # The terms of the charge-resolved sum pair as ρ and 1/ρ with equal weight when the charge
        # set is its own mirror, and 1/(1 + ρ) + 1/(1 + 1/ρ) = 1; the excitation weights pair the
        # same way, Z and Z₀ - Z being one fragmentation. The ratio of means runs both averages
        # over one set. The mean of ratios cannot: the mean of x and 1/x exceeds one.
        for label in SHIPPED_SYSTEMS
            domain = shipped_domain(label)
            A₀ = domain.system.compound.A
            @test iseven(A₀)
            @test symmetric_charge_set_is_invariant(domain) === true
            pinned = RatioCurve([A₀ ÷ 2], [0.5], [0.0], "pinned")
            weighted = ChargeResolved(
                mean_total_excitation(TEST_MASSES, domain, synthetic_kinetic_energy(domain))
            )
            for averaging in (ChargeResolved(), weighted, RatioOfMeans())
                R_T = temperature_ratio(averaging, BSFG_MODEL, domain, pinned)
                @test only(R_T.ratio) ≈ 1 atol = 1e-12
                @test only(R_T.σ) == 0
            end
            @test only(temperature_ratio(MeanOfRatios(), BSFG_MODEL, domain, pinned).ratio) < 1
        end
    end

    @testset "a single charge per mass gives every averaging the closed form" begin
        # With one fragmentation per mass the charge distribution has nothing to reduce, and all
        # three orders collapse onto R_T = [(1 - r_ν)/(ρ r_ν)]^(1/2), eq. (4) of the paper.
        domain = shipped_domain("Cf252_sf"; charges_per_mass = 1)
        @test length(domain) == length(domain.heavy_masses)
        r_ν = RatioCurve(
            collect(domain.heavy_masses),
            reference_ratio.(domain.heavy_masses),
            fill(missing, length(domain.heavy_masses)),
            "reference",
        )
        closed = Dict{Int,Float64}()
        for entry in domain
            ρ =
                level_density_parameter(BSFG_MODEL, entry.light) /
                level_density_parameter(BSFG_MODEL, entry.heavy)
            r = reference_ratio(entry.heavy.A)
            closed[entry.heavy.A] = sqrt((1 - r) / (ρ * r))
        end
        for averaging in (ChargeResolved(), RatioOfMeans(), MeanOfRatios())
            R_T = temperature_ratio(averaging, BSFG_MODEL, domain, r_ν)
            @test R_T.A_H == r_ν.A_H
            @test all(
                isapprox(R, closed[A]; rtol = 1e-12) for (A, R) in zip(R_T.A_H, R_T.ratio)
            )
        end
    end

    @testset "uncertainties propagate with the exact slope of the inverse" begin
        domain = shipped_domain("Cf252_sf")
        A_H, r, σ_r = 132, 0.3, 0.02
        curve = RatioCurve([A_H], [r], [σ_r], "synthetic")
        for averaging in (ChargeResolved(), RatioOfMeans())
            result = temperature_ratio(averaging, BSFG_MODEL, domain, curve)
            # A central difference of the inverse itself.
            h = 1e-6
            slope =
                (
                    temperature_ratio(averaging, BSFG_MODEL, domain, A_H, r + h) -
                    temperature_ratio(averaging, BSFG_MODEL, domain, A_H, r - h)
                ) / (2h)
            @test only(result.σ) ≈ abs(slope) * σ_r rtol = 1e-6
        end
        # For an effective ratio the slope is the closed form −1/(2 R_T R_a r²) of the paper.
        R_a = level_density_ratio(RatioOfMeans(), BSFG_MODEL, domain)[A_H]
        result = temperature_ratio(RatioOfMeans(), BSFG_MODEL, domain, curve)
        @test only(result.σ) ≈ σ_r / (2 * only(result.ratio) * R_a * r^2) rtol = 1e-12
    end

    @testset "singular, unquoted and uncovered points" begin
        domain = shipped_domain("Cf252_sf")
        curve = RatioCurve(
            [120, 130, 132, 140],
            [0.4, 1.0, 0.4, 0.45],
            [0.01, 0.01, missing, 0.01],
            "synthetic",
        )
        result = temperature_ratio(ChargeResolved(), BSFG_MODEL, domain, curve)
        # 120 lies below the domain and r_ν = 1 is the singular point of the relation: both
        # omitted, never extrapolated or clipped.
        @test result.A_H == [132, 140]
        # An unquoted uncertainty stays unquoted.
        @test ismissing(result.σ[1])
        @test result.σ[2] > 0
    end

    @testset "an extracted curve carries the covariance of its temperature ratio" begin
        domain = shipped_domain("Cf252_sf")
        A_H = collect(126:150)
        value = reference_ratio.(A_H) .+ 0.002 .* iseven.(A_H)
        fit = fit_segments(
            A_H, value, fill(0.004, length(A_H)); max_segments = 3, pinned_value = 0.5
        )
        data = RatioCurve(A_H, value, fill(0.004, length(A_H)), "synthetic")
        curve = ExtractedCurve("synthetic", fit, data, ChargeResolved(), BSFG_MODEL, domain)

        @test curve.coverage ≈ length(A_H) / length(domain.heavy_masses)
        @test curve.R_T.A_H == A_H
        @test size(curve.R_T_covariance) == (length(A_H), length(A_H))
        # The tabulated pointwise uncertainty is the diagonal of the propagated covariance.
        for index in eachindex(A_H)
            @test sqrt(curve.R_T_covariance[index, index]) ≈ curve.R_T.σ[index] atol = 1e-12
        end
        # Exact at the pin: zero uncertainty there, and R_T = 1 by the identity.
        @test curve.R_T.σ[1] == 0
        @test curve.R_T.ratio[1] ≈ 1 atol = 1e-12

        yields = MassYield(
            A_H, exp.(-((A_H .- 140) ./ 5) .^ 2), fill(missing, length(A_H)), "y", ""
        )
        value_cov, σ_cov = total_average(curve, yields)
        value_ind, σ_ind = total_average(curve.R_T, yields)
        # Same value; the covariance form carries the correlation the independent form drops.
        @test value_cov ≈ value_ind
        @test σ_cov > σ_ind

        mean, σ = range_mean(curve)
        @test mean ≈ sum(curve.R_T.ratio) / length(A_H)
        @test σ > 0
    end

    @testset "the jackknife covariance is the Gram matrix of the scaled deviations" begin
        refits = [
            RatioCurve([130, 131, 132, 133], [1.0, 2.0, 4.0, 7.0], fill(0.1, 4), "a"),
            RatioCurve([130, 131, 132], [2.0, 2.0, 1.0], fill(0.1, 3), "b"),
            RatioCurve([130, 131], [3.0, 5.0], fill(0.1, 2), "c"),
        ]
        Σ = FissionTemperatureRatio._jackknife_covariance([130, 131, 132, 133], refits)
        # Worked by hand. At 130 the values 1, 2, 3 deviate from their mean by −1, 0, 1, and at
        # 131 the values 2, 2, 5 by −1, −1, 2, each scaled by √(2/3); at 132 two refits hold 4
        # and 1, which deviate by ±1.5, scaled by √(1/2); at 133 one refit alone holds a value.
        @test Σ ≈ [
            4 / 3 2 -sqrt(3) / 2 0
            2 4 0 0
            -sqrt(3) / 2 0 9 / 4 0
            0 0 0 0
        ] atol = 1e-14
        @test Σ == Σ'
        @test FissionTemperatureRatio._jackknife_covariance([130, 131], RatioCurve[]) ==
            zeros(2, 2)
    end

    @testset "a curve with refits carries the larger of its propagated uncertainty and the jackknife" begin
        domain = shipped_domain("Cf252_sf")
        A_H = collect(126:150)
        base = reference_ratio.(A_H) .+ 0.002 .* iseven.(A_H)
        function extracted(value, label)
            fit = fit_segments(
                A_H, value, fill(0.004, length(A_H)); max_segments = 3, pinned_value = 0.5
            )
            data = RatioCurve(A_H, value, fill(0.004, length(A_H)), label)
            return ExtractedCurve(label, fit, data, ChargeResolved(), BSFG_MODEL, domain)
        end
        # The refits leave the pin where it is and differ from one another by a slope, far more
        # than the propagated uncertainty of the curve.
        tilted(shift) = extracted(base .+ shift .* (A_H .- 126) ./ 24, "refit")
        curve = extracted(base, "synthetic")
        refits = [tilted(shift) for shift in (-0.05, 0.0, 0.08)]
        @test all(refit.R_T.A_H == curve.R_T.A_H for refit in refits)
        @test all(refit.r_ν.A_H == curve.r_ν.A_H for refit in refits)
        raised = ExtractedCurve(curve, refits)

        @test curve.R_T_jackknife_covariance === nothing
        @test raised.label == curve.label
        @test raised.fit === curve.fit
        @test raised.coverage == curve.coverage
        @test raised.r_ν.ratio == curve.r_ν.ratio
        @test raised.R_T.ratio == curve.R_T.ratio
        @test raised.R_T_covariance == curve.R_T_covariance
        # The delete-one jackknife standard error of k values.
        jackknife(θ) = sqrt((length(θ) - 1) / length(θ) * sum(abs2, θ .- sum(θ) / length(θ)))
        for index in eachindex(A_H)
            σ_r = jackknife([refit.r_ν.ratio[index] for refit in refits])
            σ_R = jackknife([refit.R_T.ratio[index] for refit in refits])
            @test raised.r_ν.σ[index] ≈ max(curve.r_ν.σ[index], σ_r) atol = 1e-14
            @test raised.R_T.σ[index] ≈ max(curve.R_T.σ[index], σ_R) atol = 1e-14
            @test sqrt(raised.R_T_jackknife_covariance[index, index]) ≈ σ_R atol = 1e-14
        end
        # Every refit holds the pin: no spread there, and the uncertainty stays zero.
        @test raised.R_T.σ[1] == 0
        @test raised.R_T.σ[end] > 2 * curve.R_T.σ[end]

        # An average of the curve takes the jackknife of the same average over the refits.
        yields = MassYield(
            A_H, exp.(-((A_H .- 140) ./ 5) .^ 2), fill(missing, length(A_H)), "y", ""
        )
        value, σ = total_average(raised, yields)
        σ_jackknife = jackknife([first(total_average(refit, yields)) for refit in refits])
        @test value == first(total_average(curve, yields))
        @test σ_jackknife > last(total_average(curve, yields))
        @test σ ≈ σ_jackknife rtol = 1e-10
        σ_mean = last(range_mean(raised))
        @test σ_mean > last(range_mean(curve))
        @test σ_mean ≈ jackknife([first(range_mean(refit)) for refit in refits]) rtol = 1e-10

        # Where the propagated uncertainty is the larger it stands: refits that barely differ.
        unraised = ExtractedCurve(curve, [tilted(shift) for shift in (-1e-6, 0.0, 1e-6)])
        @test unraised.R_T.σ == curve.R_T.σ
        @test unraised.r_ν.σ == curve.r_ν.σ
        @test last(total_average(unraised, yields)) == last(total_average(curve, yields))

        @test_throws ArgumentError ExtractedCurve(curve, refits[1:1])
    end
end
