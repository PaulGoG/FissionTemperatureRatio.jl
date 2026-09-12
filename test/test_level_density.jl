
@testset "level density parameter" begin
    prescription = BackShiftedFermiGas()

    @testset "magnitude is physical for fission fragments" begin
        for (A, Z) in ((100, 40), (132, 50), (140, 54))
            a = level_density_parameter(prescription, A, Z, TEST_MASSES)
            @test a isa Float64
            @test 0 < a < A / 4
        end
    end

    @testset "shell closures suppress the parameter" begin
        # The systematic carries the shell correction, so the doubly magic heavy fragment sits far
        # below the smooth trend of its neighbours. This suppression is what gives the level
        # density parameter ratio its structure around A_H = 130.
        a(A, Z) = level_density_parameter(prescription, A, Z, TEST_MASSES)
        magic = a(132, 50)
        @test magic < a(128, 50)
        @test magic < a(136, 54)
        # Reduced by a factor of about three relative to a mid-shell fragment.
        @test magic / 132 < 0.5 * a(100, 40) / 100
    end

    @testset "absent mass excesses propagate as missing" begin
        @test ismissing(level_density_parameter(prescription, 400, 150, TEST_MASSES))
        @test ismissing(shell_correction(400, 150, TEST_MASSES))
    end

    @testset "scaling with mass number" begin
        # At fixed shell correction the systematic is a power law in A, so a heavier nuclide has
        # the larger parameter.
        light = level_density_parameter(prescription, 100, 40, TEST_MASSES)
        heavy = level_density_parameter(prescription, 150, 60, TEST_MASSES)
        @test heavy > light
    end
end
