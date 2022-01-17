using DrillMudsThermalProps
using Test

@testset "DrillMudsThermalProps.jl" begin
    # Write your tests here.

    water=Water()

    @test water.ρ==1000.0
    


end
