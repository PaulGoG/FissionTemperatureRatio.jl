
DATA_AVAILABLE && @testset "level density parameter" begin
    prescription = BSFG_PRESCRIPTION

    @testset "magnitude is physical for fission fragments" begin
        for (A, Z) in ((100, 40), (132, 50), (140, 54))
            a = level_density_parameter(prescription, A, Z)
            @test a isa Float64
            @test 0 < a < A / 4
        end
    end

    @testset "shell closures suppress the parameter" begin
        # The systematic carries the shell correction, so the doubly magic heavy fragment sits far
        # below the smooth trend of its neighbours. This suppression is what gives the level
        # density parameter ratio its structure around A_H = 130.
        a(A, Z) = level_density_parameter(prescription, A, Z)
        magic = a(132, 50)
        @test magic < a(128, 50)
        @test magic < a(136, 54)
        # Reduced by a factor of about three relative to a mid-shell fragment.
        @test magic / 132 < 0.5 * a(100, 40) / 100
    end

    @testset "absent mass excesses propagate as missing" begin
        @test ismissing(level_density_parameter(prescription, 400, 150))
        @test ismissing(shell_correction(400, 150, TEST_MASSES))
    end

    @testset "scaling with mass number" begin
        # At fixed shell correction the systematic is a power law in A, so a heavier nuclide has
        # the larger parameter.
        light = level_density_parameter(prescription, 100, 40)
        heavy = level_density_parameter(prescription, 150, 60)
        @test heavy > light
    end

    @testset "Gilbert-Cameron" begin
        # Provided for assessing the dependence on the prescription. Away from closed shells it
        # returns parameters well above the back-shifted Fermi gas, which is the documented reason
        # it is not the default. At the doubly magic heavy fragment both are suppressed, so the
        # comparison is made where the two prescriptions actually part company.
        for (A, Z) in ((100, 40), (140, 54), (150, 60))
            gc = level_density_parameter(GC_PRESCRIPTION, A, Z)
            bsfg = level_density_parameter(BSFG_PRESCRIPTION, A, Z)
            @test gc isa Float64
            @test gc > bsfg
        end
        @test level_density_parameter(GC_PRESCRIPTION, 132, 50) <
            level_density_parameter(GC_PRESCRIPTION, 140, 54)
        # Outside the tabulated nucleon numbers the prescription is undefined.
        @test ismissing(level_density_parameter(GC_PRESCRIPTION, 400, 200))
    end
end
