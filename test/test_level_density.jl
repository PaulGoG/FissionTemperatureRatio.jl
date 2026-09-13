
DATA_AVAILABLE && @testset "level density parameter" begin
    model = BSFG_MODEL

    @testset "magnitude is physical for fission fragments" begin
        for (A, Z) in ((100, 40), (132, 50), (140, 54))
            a = level_density_parameter(model, A, Z)
            @test a isa Float64
            @test 0 < a < A / 4
        end
    end

    @testset "shell closures suppress the parameter" begin
        # The systematic carries the shell correction, so the doubly magic heavy fragment sits far
        # below the smooth trend of its neighbours. This suppression is what gives the level
        # density parameter ratio its structure around A_H = 130.
        a(A, Z) = level_density_parameter(model, A, Z)
        magic = a(132, 50)
        @test magic < a(128, 50)
        @test magic < a(136, 54)
        # Reduced by a factor of about three relative to a mid-shell fragment.
        @test magic / 132 < 0.5 * a(100, 40) / 100
    end

    @testset "absent mass excesses propagate as missing" begin
        @test ismissing(level_density_parameter(model, 400, 150))
        @test ismissing(shell_correction(400, 150, TEST_MASSES))
    end

    @testset "scaling with mass number" begin
        # At fixed shell correction the systematic is a power law in A, so a heavier nuclide has
        # the larger parameter.
        light = level_density_parameter(model, 100, 40)
        heavy = level_density_parameter(model, 150, 60)
        @test heavy > light
    end

    @testset "Gilbert-Cameron" begin
        # Provided for assessing the dependence on the model. Away from closed shells it
        # returns parameters well above the back-shifted Fermi gas, which is the documented reason
        # it is not the default. At the doubly magic heavy fragment both are suppressed, so the
        # comparison is made where the two models actually part company.
        for (A, Z) in ((100, 40), (140, 54), (150, 60))
            gc = level_density_parameter(GC_MODEL, A, Z)
            bsfg = level_density_parameter(BSFG_MODEL, A, Z)
            @test gc isa Float64
            @test gc > bsfg
        end
        @test level_density_parameter(GC_MODEL, 132, 50) <
            level_density_parameter(GC_MODEL, 140, 54)
        # The table declares `n S(N) S(Z)`, and the expression looks S(Z) up at the proton number
        # and S(N) at the neutron number. Exchanging the two columns is not absorbed by their sum:
        # it evaluates the neutron correction at the proton number and the reverse, silently. The
        # first tabulated row fixes the order, and the closed form fixes the lookup keys.
        @test TEST_SHELLS.S_N[11] ≈ 6.80
        @test TEST_SHELLS.S_Z[11] ≈ -2.91
        # Against Table III of Gilbert and Cameron, Can. J. Phys. 43, 1446 (1965), at the nucleon
        # numbers that matter most for fission fragments: the proton and neutron shell closures,
        # and the heaviest entries the table carries.
        for (n, S_Z, S_N) in (
            (50, -19.83, 12.88),
            (52, -18.35, 13.71),
            (54, -16.54, 15.16),
            (82, -8.86, 9.09),
            (98, -7.74, 9.65),
        )
            @test TEST_SHELLS.S_Z[n] ≈ S_Z
            @test TEST_SHELLS.S_N[n] ≈ S_N
        end
        # Eq. (20), the undeformed correlation, not Eq. (21) which offsets by 0.120 instead.
        @test level_density_parameter(GC_MODEL, 132, 50) ≈
            132 * (0.00917 * (TEST_SHELLS.S_Z[50] + TEST_SHELLS.S_N[82]) + 0.142)
        @test level_density_parameter(GC_MODEL, 132, 50) ≈
            132 * (9.17e-3 * (TEST_SHELLS.S_Z[50] + TEST_SHELLS.S_N[82]) + 1.42e-1)
        # Outside the tabulated nucleon numbers the model is undefined.
        @test ismissing(level_density_parameter(GC_MODEL, 400, 200))
    end
end
