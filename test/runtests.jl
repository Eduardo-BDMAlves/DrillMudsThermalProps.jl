using DrillMudsThermalProps
using Test

@testset "Entire DrillMudsThermalProps.jl" begin
    atol_tests=1.0E-7
    # Write your tests here.
    include("WaterTests.jl")
    include("Solids.jl")
    include("MixtureDensity.jl")



    @test isa(try 
                rho(101325.0,298.15,TestFluid()) 
            catch e
                e 
            end,Exception)

end