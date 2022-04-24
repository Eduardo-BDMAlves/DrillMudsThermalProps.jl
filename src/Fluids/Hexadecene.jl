export Hexadecene,
        EOS,
        PRSV,
        PR,
        RK


abstract type EOS end


struct PRSV end
struct PR end
struct RK end



struct Hexadecene <: Fluid


    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    MM::Float64

    Tc::Float64
    Pc::Float64
    ω::Float64
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
            722.0124308882724,
            1370.0E3,
            0.693
            # 730.0
        )

    end

    #https://wtt-pro.nist.gov/wtt-pro/index.html?cmp=1-hexadecene#1-hexadecene;06;0H;0A;4g;1g;7g;3g;5g;Hg;8g,298.15,413.15,10;8b;O6;OA;OH;Qg;TK;TC;T2;Tg;Re;zK;zC;z2;zg;Ye;aK;aC;a2;ag;ae;hf;mg;jg;ja;na;ng;ya;yf/c;0,0/a;782,144/8gR;805,301/8gT;72,426,508,382/

end



function rho(P,T,fluid::Hexadecene,EOS::PRSV)
    # return fluid.ρ
    P=P/1000
    ## Based in the PRSV EOS
    k0=1.3850916#/0.8
    k1=-0.048310908#/1000

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

    # f(v) = v^3+A*v^2+B*v+C

    # v=find_zero(f,(fluid.MM/fluid.ρ/2,fluid.MM/fluid.ρ*2),A42(),xtol=1.0E-8)

    C1=A
    C2=B
    C3=C
    C4=C1/3

    D=0.5*(C2*C4-C3)-C4^3
    E=C4^2-C2/3
    Δ=D^2-E^3

    tol=1.0E-8
    v=0.0
    # println((
    #     Δ=Δ,
    #     E=E,
    #     D=D,
    #     C4=C4,
    #     C3=C3,
    #     C2=C2,
    #     C1=C1
    # ))

    if Δ≥tol
        F=cbrt(D+sqrt(Δ))
        G=cbrt(D-sqrt(Δ))
        # H=-(0.5*(F+G)-C4)
        v=F+G-C4

    elseif Δ≥-tol
        v1=2cbrt(D)-C4
        v2=-cbrt(D)-C4
        # println((v1=v1,v2=v2))
        if 0.23 ≤ v1 ≤ 0.38
            v=v1
        elseif 0.23 ≤ v2 ≤ 0.38
            v=v2
        else
            v=0.0
        end
    else
        J=acos(D/sqrt(E^3))
        v1=2sqrt(E)*cos(J/3)-C4
        v2=2sqrt(E)*cos(J/3 + 2π/3)-C4
        v3=2sqrt(E)*cos(J/3+4π/3)-C4

        # println((v1=v1,v2=v2,v3=v3))

        if 0.23 ≤ v1 ≤ 0.38
            v=v1
        elseif 0.23 ≤ v2 ≤ 0.38
            v=v2
        elseif 0.23 ≤ v3 ≤ 0.38
            v=v3
        else
            v=0.0
        end
    end

    # println((v=v,))


    # v=find_zeros(f,(fluid.MM/fluid.ρ/10,fluid.MM/fluid.ρ*10))
    # return fluid.MM/v


    Av=0.808
    Bv=2.575

    Dv=-0.04279848596094332#-0.6144070066691208*1000#brute force...

    v_new=(v+Dv)*(1-exp(Av-Bv/Tr^3))
    # v_new=v
    return fluid.MM/v_new
end


function rho(P,T,fluid::Hexadecene,EOS::RK)

    R=8.31446261815324

    a=0.42748*R*R*fluid.Tc^(5/2)/fluid.Pc

    b=0.08664*R*fluid.Tc/fluid.Pc

    Tr=T/fluid.Tc

    alpha=(1+(0.37464+1.54226fluid.ω-0.26992fluid.ω^2)*(1-sqrt(Tr)))^2

    println((α=alpha,a=a,b=b))



    

    # A=a*alpha*P/(R*T)^2
    # B=b*P/(R*T)


    f(v)=R*T/(v-b)-a/(sqrt(T)*v*(v+b))
    # C1=(B-1)
    # C2=A-2B-3B^2
    # C3=(B^2+B-A)*B

    # f(v) = v^3+C1*v^2+C2*v+C3

    println((a_eq=fluid.MM/fluid.ρ/10,b_eq=fluid.MM/fluid.ρ*4))
    # v=find_zero(f,(fluid.MM/fluid.ρ/10000000,fluid.MM/fluid.ρ),A42(),xtol=1.0E-8)
    v=find_zero(f,fluid.MM/fluid.ρ,tol=1.0E-8)

    return fluid.MM/v
end


function rho(P,T,fluid::Hexadecene,EOS::PR)

    R=8.31446261815324

    a=0.45724(R*fluid.Tc)^2/fluid.Pc

    b=0.07780R*fluid.Tc/fluid.Pc

    Tr=T/fluid.Tc

    alpha=(1+(0.37464+1.54226fluid.ω-0.26992fluid.ω^2)*(1-sqrt(Tr)))^2

    println((α=alpha,a=a,b=b))

    A=a*alpha*P/(R*T)^2
    B=b*P/(R*T)

    C1=(B-1)
    C2=A-2B-3B^2
    C3=(B^3+B^2-A*B)

    C4=C1/3

    D=0.5*(C2*C4-C3)-C4^3
    E=C4^2-C2/3
    Δ=D^2-E^3

    tol=1.0E-6
    v=0.0
    println((
        Δ=Δ,
        E=E,
        D=D,
        C4=C4,
        C3=C3,
        C2=C2,
        C1=C1
    ))

    if Δ≥tol
        F=cbrt(D+sqrt(Δ))
        G=cbrt(D-sqrt(Δ))
        # H=-(0.5*(F+G)-C4)
        v=F+G-C4

    elseif Δ≥-tol
        v1=2cbrt(D)-C4
        v2=-cbrt(D)-C4
        println((v1=v1,v2=v2))
        if 0.23 ≤ v1 ≤ 0.38
            v=v1
        elseif 0.23 ≤ v2 ≤ 0.38
            v=v2
        else
            v=0.0
        end
    else
        J=acos(D/sqrt(E^3))
        v1=2sqrt(E)*cos(J/3)-C4
        v2=2sqrt(E)*cos(J/3 + 2π/3)-C4
        v3=2sqrt(E)*cos(J/3+4π/3)-C4

        println((v1=v1,v2=v2,v3=v3))

        if 0.23 ≤ v1 ≤ 0.38
            v=v1
        elseif 0.23 ≤ v2 ≤ 0.38
            v=v2
        elseif 0.23 ≤ v3 ≤ 0.38
            v=v3
        else
            v=0.0
        end
    end
    
    v_m=v*R*T/P

    return fluid.MM/v_m
    return fluid.MM/v


    # return 0.5
end

rho(P,T,fluid::Hexadecene)=rho(P,T,fluid,PR())




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






