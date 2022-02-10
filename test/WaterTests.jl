@testset "Water tests." begin

    water = Water()


    Ts = collect([273.153:10.0:140.3+273.15]...)
    Ps = collect([[101325.0]..., 5.0E6:19.0E6:100.0E6...])


    rhos_CP = [
        999.843288668736
        999.7022064047571
        998.206531114165
        995.648547727687
        992.2152054029182
        988.0336893152614
        983.1942808217368
        977.7629140511095
        971.7885282590977
        965.3075721956923
        1002.321301004123
        1002.0310410719416
        1000.4389567920352
        997.8207863375475
        994.3540212626696
        990.1600120842329
        985.325295236465
        979.9133791196485
        973.9716718278099
        967.5357377570548
        960.6319955267068
        953.2794593328697
        945.4908642455383
        937.2733731926065
        928.628980933184
    ]

    Cps_CP = [
        4219.434432031926
        4195.15402061601
        4184.048832139037
        4179.819080630115
        4179.415075226967
        4181.343152236771
        4184.954590522909
        4190.068861555105
        4196.75552379441
        4205.208431029631
        4195.6772230196875
        4176.951566752204
        4169.021633960576
        4166.726174432052
        4167.529584370848
        4170.198442209381
        4174.226023012843
        4179.512509552163
        4186.1735172851095
        4194.427743361589
        4204.535630211352
        4216.770290859532
        4231.408130310603
        4248.731120668077
        4269.035936461646
    ]

    ks_CP = [
        0.5556573479629365
        0.57878368629924
        0.598017660579779
        0.6143967504241387
        0.6284896176744261
        0.6406244513552712
        0.6510031477223245
        0.6597606493966989
        0.6669962636640698
        0.6727901194178599
        0.5593479468610951
        0.5819555952701735
        0.6008790707937224
        0.6170777712382046
        0.6310751274302673
        0.6431717005931662
        0.653551885098157
        0.6623392583494603
        0.6696256425276464
        0.6754862275416963
        0.6799871759111179
        0.6831891713047554
        0.6851488181656393
        0.6859189511453694
        0.6855484362606141
    ]

    mus_CP=[
    0.0017915688995613212
    0.0013057868784517098
    0.0010015225439034675
    0.0007971708792954054
    0.0006526919175727183
    0.0005464887396638285
    0.00046601393786597313
    0.0004035315809411959
    0.0003540373865383654
    0.00031416450875558944
    0.0017808395544338116
    0.0013014716443436355
    0.0010000705264838394
    0.0007971003389448608
    0.0006533214417919576
    0.0005474817473247409
    0.00046719513565785306
    0.00040480615876922705
    0.00035535285387594236
    0.0003154916084978094
    0.0002828958170615312
    0.0002559024852849985
    0.00023329557520111158
    0.00021416827628534935
    0.0001978326774638032
    ]

    rhos = Float64[]
    Cps = Float64[]
    ks = Float64[]
    mus = Float64[]

    ips = [
        (101325.0, 273.153)
        (101325.0, 283.153)
        (101325.0, 293.153)
        (101325.0, 303.153)
        (101325.0, 313.153)
        (101325.0, 323.153)
        (101325.0, 333.153)
        (101325.0, 343.153)
        (101325.0, 353.153)
        (101325.0, 363.153)
        (5.0e6, 273.153)
        (5.0e6, 283.153)
        (5.0e6, 293.153)
        (5.0e6, 303.153)
        (5.0e6, 313.153)
        (5.0e6, 323.153)
        (5.0e6, 333.153)
        (5.0e6, 343.153)
        (5.0e6, 353.153)
        (5.0e6, 363.153)
        (5.0e6, 373.153)
        (5.0e6, 383.153)
        (5.0e6, 393.153)
        (5.0e6, 403.153)
        (5.0e6, 413.153)
    ]

    for (P,T) ∈ ips
            rho = DrillMudsThermalProps.rho(P, T, water)
            cp = DrillMudsThermalProps.Cp(P, T, water)
            k = DrillMudsThermalProps.therm_cond(P, T, water)
            mu = DrillMudsThermalProps.visc(P,T,water)

            append!(rhos, rho)
            append!(Cps, cp)
            append!(ks, k)
            append!(mus,mu)
    end


    dev_rho = abs.(rhos_CP .- rhos) ./ rhos_CP

    @testset "Density" begin
        for (dev,ip) in zip(dev_rho,ips)
            @test (dev ≤ 0.01)
        end
    end

    dev_Cp=abs.(Cps_CP.-Cps)./Cps
    @testset "Cp" begin
        for (dev,ip) in zip(dev_Cp,ips)
            @test (dev ≤ 0.05)
        end
    end

    dev_mu=abs.(mus.-mus_CP)./mus_CP
    @testset "Viscosity" begin
        for (dev,ip) in zip(dev_mu,ips)
            @test (dev ≤ 0.07)
        end
    end



    dev_ks=abs.(ks-ks_CP)./ks_CP
    @testset "therm. cond." begin
        for (dev,ip) in zip(dev_ks,ips)
            @test (dev ≤ 0.05)
        end
    end

end





















