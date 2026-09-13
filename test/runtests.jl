using Test
using FissionTemperatureRatio

include("fixtures.jl")

@testset "FissionTemperatureRatio" begin
    include("test_quality.jl")
    include("test_mass_data.jl")
    include("test_level_density.jl")
    include("test_fragmentation.jl")
    include("test_multiplicity_ratio.jl")
    include("test_yields.jl")
    include("test_temperature_ratio.jl")
    include("test_segmented_fit.jl")
    include("test_configuration.jl")
end
