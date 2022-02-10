@testset "Mixture density" begin

    @testset "Limit concentrations leads to pure substance?" begin
        water = Water()
        brine860 = BrineNaCl(8.6)
        brine920 = BrineNaCl(9.2)
        barite = Barite()

        Ps=LinRange(101325.0,100.0E6,5)
        Ts=LinRange(283.15,413.15,5)
        @testset "Water + Barite mixture" begin
            FList=[water]
            SLits=[barite]
            fs=[1.0,0.0]

            WatBar0=DrillFluid(FList,SLits,fs=fs)

            fs=[0.0,1.0]
            WatBar100=DrillFluid(FList,SLits,fs=fs)

            for P in Ps
                for T in Ts
                    rho_water = rho(P, T, water)
                    rho_barite = rho(P, T, barite)
                    rho_mix0 = rho(P, T, WatBar0)
                    rho_mix100 = rho(P, T, WatBar100)
                
                    @test rho_water == rho_mix0
                
                    @test rho_barite == rho_mix100
                
                end
            end
            
        end

        @testset "Brine 8.6ppg + Barite mixture" begin
            FList = [brine860]
            SLits = [barite]
            fs = [1.0, 0.0]

            BrBar0 = DrillFluid(FList, SLits, fs = fs)

            fs = [0.0, 1.0]
            BrBar100 = DrillFluid(FList, SLits, fs = fs)

            for P in Ps
                for T in Ts
                    rho_Br = rho(P, T, brine860)
                    rho_barite = rho(P, T, barite)
                    rho_mix0 = rho(P, T, BrBar0)
                    rho_mix100 = rho(P, T, BrBar100)

                    @test rho_Br == rho_mix0

                    @test rho_barite == rho_mix100

                end
            end

        end


        @testset "Brine 9.2ppg + Barite mixture" begin
            FList = [brine920]
            SLits = [barite]
            fs = [1.0, 0.0]

            BrBar0 = DrillFluid(FList, SLits, fs = fs)

            fs = [0.0, 1.0]
            BrBar100 = DrillFluid(FList, SLits, fs = fs)

            for P in Ps
                for T in Ts
                    rho_Br = rho(P, T, brine920)
                    rho_barite = rho(P, T, barite)
                    rho_mix0 = rho(P, T, BrBar0)
                    rho_mix100 = rho(P, T, BrBar100)

                    @test rho_Br == rho_mix0

                    @test rho_barite == rho_mix100

                end
            end

        end

    end
    
end