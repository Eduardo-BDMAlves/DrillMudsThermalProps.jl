export Hexadecene

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



function rho(P,T,fluid::Hexadecene)
    return fluid.ρ
end

function Cp(P,T,fluid::Hexadecene)
    return fluid.cₚ
end









