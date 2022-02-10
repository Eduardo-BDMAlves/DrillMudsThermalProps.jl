@testset "Solids tests" begin
    @testset "Barite" begin
        barite=Barite()        

        @test barite.ρ≈rho(101325.0,298.15,barite)
        @test barite.cₚ≈Cp(101325.0,298.15,barite)
        @test barite.k≈therm_cond(101325.0,298.15,barite)
    end

    @testset "Hematite" begin
        hematite=Hematite()        

        @test hematite.ρ≈rho(101325.0,298.15,hematite)
        @test hematite.cₚ≈Cp(101325.0,298.15,hematite)
        @test hematite.k≈therm_cond(101325.0,298.15,hematite)
    end

    @testset "CaCO3" begin
        carbonite=CaCO3()

        @test carbonite.ρ≈rho(101325.0,298.15,carbonite)
        @test carbonite.cₚ≈Cp(101325.0,298.15,carbonite)
        @test carbonite.k≈therm_cond(101325.0,298.15,carbonite)
    end


end