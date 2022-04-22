export Hexadecene

struct Hexadecene <: Fluid


    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    MM::Float64

    Tc::Float64
    #viscosity ver no arquivo c16.pdf

    function Hexadecene()
        MM=224.43

        new(
            :Hexadecene,
            781.1,#no Nist 777.91
            # 679.9173943031689,
            357.2/MM*1.0E3,
            3.0E-3,
            0.1423,
            MM,
            722.0124308882724
            # 730.0
        )

    end

    #https://wtt-pro.nist.gov/wtt-pro/index.html?cmp=1-hexadecene#1-hexadecene;06;0H;0A;4g;1g;7g;3g;5g;Hg;8g,298.15,413.15,10;8b;O6;OA;OH;Qg;TK;TC;T2;Tg;Re;zK;zC;z2;zg;Ye;aK;aC;a2;ag;ae;hf;mg;jg;ja;na;ng;ya;yf/c;0,0/a;782,144/8gR;805,301/8gT;72,426,508,382/

end



function rho(P,T,fluid::Hexadecene)
    # return fluid.ρ

    ## Based in the PRSV EOS
    k0=1.3850916
    k1=-0.048310908

    ac=11241.446#*1.0E6
    b=0.31861633#*1.0E3

    Tr=T/fluid.Tc
    # Tr=T/(0.020464*ac/b)
    κ=k0+k1*(1+sqrt(Tr))*(0.7-Tr)
    alpha=(1+κ*(1-sqrt(Tr)))*(1+κ*(1-sqrt(Tr)))

    R=8.31446261815324

    A=b-R*T/P
    B=ac*alpha/P-2R*T*b/P-3b^2
    C=(b^2+R*T*b/P-ac*alpha/P)*b

    f(v) = v^3+A*v^2+B*v+C

    v=find_zero(f,(fluid.MM/fluid.ρ/2,fluid.MM/fluid.ρ*2),A42(),xtol=1.0E-8)

    # v=find_zeros(f,(fluid.MM/fluid.ρ/10,fluid.MM/fluid.ρ*10))
    return fluid.MM/v


    # Av=0.808
    # Bv=2.575

    # Dv=-0.6144070066691208*1000#brute force...

    # v_new=(v+Dv)*(1-exp(Av-Bv/Tr^3))

    # return fluid.MM/v_new
end



function visc(P,T,fluid::Hexadecene)
    @warn "Not implemented, using constant value for implementation of mixture... Implement the rule."
    return fluid.μ
end

function Cp(P,T,fluid::Hexadecene)
    return fluid.cₚ
end


function therm_cond(P,T,fluid::Hexadecene)
    # @warn "Not implemented, using constant value for implementation of mixture... Implement the rule."
    # return fluid.k
    A=-1.7204
    B=1.0172
    C=722.0
    log_lamb=A+B*(1.0-T/C)^(2/7)
    return 10^log_lamb
end






