using DrillMudsThermalProps
using Test

@testset "Entire DrillMudsThermalProps.jl" begin
    # Write your tests here.
    include("WaterTests.jl")
    include("Solids.jl")
    include("MixtureDensity.jl")

end