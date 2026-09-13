
@testset "mass excess table" begin
    if DATA_AVAILABLE
        @test mass_excess(TEST_MASSES, 132, 50) isa Float64
        @test ismissing(mass_excess(TEST_MASSES, 400, 50))
    end

    @testset "reader rejects malformed input" begin
        @test_throws ArgumentError read_mass_excess_table(joinpath(@__DIR__, "absent.dat"))
        mktempdir() do directory
            path = joinpath(directory, "no_proton.dat")
            write(path, "50 132 Sn -76.5 0.3\n")
            @test_throws ArgumentError read_mass_excess_table(path)
        end
    end
end
