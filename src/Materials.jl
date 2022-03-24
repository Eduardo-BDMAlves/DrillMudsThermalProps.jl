export TestSolid,
    TestFluid,
    Water,
    Barite,
    Fluid,
    Solid,
    Fluid_sets,
    BrineNaCl,
    CaCO3,
    Hematite,
    Hexadecene


abstract type Material end

abstract type Solid <: Material end
abstract type Fluid <: Material end

# Fluid_sets=Union{(Fluid, (Fluid)...)...}

# struct Water <: Fluid end

struct TestSolid <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    function TestSolid()
        new(
            :TestSolid,
            1000.0,
            4000.0,
            1.002E-3,
            0.58
        )
    end
end


struct TestFluid <: Fluid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    function TestFluid()
        new(
            :Water,
            1000.0,
            4000.0,
            1.002E-3,
            0.58
        )
    end
end




struct Water <: Fluid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    function Water()
        new(
            :Water,
            1000.0,
            4000.0,
            1.002E-3,
            0.58
        )
    end
end


struct BrineNaCl <: Fluid

    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64
    x::Float64



    function BrineNaCl(ppg = 9.8)

        ppgs = 8.6:0.1:9.8
        xs = [0.45501, 0.061330, 0.077159, 0.092986, 0.108814, 0.124642, 0.140470, 0.156297, 0.172125, 0.187953, 0.203782, 0.219610, 0.235437]

        itp = LinearInterpolation(ppgs, xs)

        x = itp(ppg)

        P = 101325.0
        T = 273.15 + 20.0

        #rho
        a, b, c, d, g, h, m, n, p, q, r, s = 8.23073018e+02, 1.38381011e+00, -2.85110020e-03, -9.11569454e-18, 4.04861229e-08, 2.38570427e-15, -2.65736884e-09, 9.52410475e-12, 1.27087646e+03, -7.83633987e-01, 6.26293103e-09, -2.11977575e-03

        d1 = a + b * T + c * T^2
        d2 = g * P + h * P^2
        d3 = m * P * T + n * P * T^2 + d * P^2 * T
        d4 = p * x + q * x * T + r * x * P * T + s * x * T^2

        rho_std = d1 + d2 + d3 + d4

        # Cp
        b, c, d, e, f, g, i, j, k, l, m = (2.48218846e+03, 3.80939153e-06, -3.45541918e-15, -2.57191794e+03, 3.85762099e+03, 1.93733545e-07, 2.27236796e+00, -9.07624069e-09, 1.27273598e-17, 1.38313250e+00, -1.96832796e+00)

        h_1 = b + c * P + d * P^2 + e * x + f * x^2 + g * x * P
        h_2 = 2 * (i + j * P + k * P^2 + l * x + m * x^2) * T
        cp_std = h_1 + h_2

        # mu 


        a, b, c, d, g, h, m, p, q, r = (-9.29821669e+00, -5.62181704e+02, -2.60785572e-07, 2.24278114e+03, -1.75402843e-16, 3.70625444e+05, -4.33694111e+05, 4.47607787e-12, -1.44053809e-20, -1.33905655e-05)


        v_1 = a + (b + c * P + d * x + g * P^2) / T
        v_2 = (h + m * x) / T^2
        v_3 = p * P * T + q * P^2 * T + r * T^2 * x^2
        mu_std = exp(v_1 + v_2 + v_3)


        # therm cond

        a, b, c, d, g, h, m, n, p, q, r, s, t = (-9.60221679e-01, 1.06170815e-02, -2.25969876e-05, 1.57190072e-08, -8.34381967e-10, -2.63526539e-19, 4.75581304e-12, -3.22756205e-15, 4.47716476e-01, -4.57135183e-03, 1.55078050e-09, -2.82355304e-12, 8.68547458e-06)

        k_1 = a + b * T + c * T^2 + d * T^3
        k_2 = g * P + h * P^2
        k_3 = m * P * T + n * P * T^2
        k_4 = p * x + q * x * T + r * x * P + s * x * P * T + t * x * T^2

        k_std = k_1 + k_2 + k_3 + k_4



        new(
            :Brine_NaCl,
            rho_std,
            cp_std,
            mu_std,
            k_std,
            x
        )
    end


end



struct Hexadecene <: Fluid


    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    MM::Float64

    #viscosity ver no arquivo c16.pdf

    function Hexadecene()
        MM=224.43

        new(
            :Hexadecene,
            781.1,#no Nist 777.91
            357.2/MM*1.0E3,
            3.0E-3,
            0.1423,
            MM
        )

    end

    #https://wtt-pro.nist.gov/wtt-pro/index.html?cmp=1-hexadecene#1-hexadecene;06;0H;0A;4g;1g;7g;3g;5g;Hg;8g,298.15,413.15,10;8b;O6;OA;OH;Qg;TK;TC;T2;Tg;Re;zK;zC;z2;zg;Ye;aK;aC;a2;ag;ae;hf;mg;jg;ja;na;ng;ya;yf/c;0,0/a;782,144/8gR;805,301/8gT;72,426,508,382/

end






## General computation common to every fluid

function (F::Fluid)(P, T)
    ρ = rho(P, T, F)
    κ = compress(P, T, F)
    β = therm_expand(P, T, F)
    C = Cp(P, T, F)
    k = therm_cond(P, T, F)
    μ = visc(P, T, F)

    println("hello")
    return (ρ = ρ,
        κ = κ,
        β = β,
        Cₚ = C,
        k = k,
        μ = μ)
end









## Solid

struct Barite <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    function Barite()
        new(
            :Barite,
            4100.0,
            4600.0,
            0.09515
        )
    end
end




struct CaCO3 <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    # http://www.matweb.com/search/datasheet_print.aspx?matguid=bea4bfa9c8bd462093d50da5eebe78ac
    # https://thermtest.com/thermal-resources/materials-database
    # Caenn -> drilling mud
    function CaCO3()
        new(
            :CaCO3,
            2610.0,
            834.0,
            2.259
        )
    end
end



struct Hematite <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C1317608&Type=JANAFS&Table=on
    # https://www.researchgate.net/publication/252121179_Thermal_Conductivity_of_Magnetite_and_Hematite
    function Hematite()
        new(
            :Hematite,
            4900.0,
            649.8627065488718,
            6.413265500000001
        )
    end
end
