export Hexadecene



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


    RK_params::NamedTuple
    PR_params::NamedTuple
    PRSV_params::NamedTuple

    TammannTait_params::NamedTuple


    #viscosity ver no arquivo c16.pdf

    function Hexadecene()
        MM=224.43

        Tc=722.0124308882724
        Pc=1370.0E3
        ω=0.767
        k0=1.3850916
        k1=-0.048310908

        RK_params=set_up_EOS(RK(),Tc,Pc,ω)
        PR_params=set_up_EOS(PR(),Tc,Pc,ω)
        PRSV_params=set_up_EOS(PRSV(),Tc,Pc,k0,k1)

        ## TammannTait_params fited from data...
        TTps=(
            E=-0.11503916234335715,
            F0=-5.0268235197563744e8,
            F1=1.6685928787241152e6,
            F2=-1468.7483988370013
        )



        new(
            :Hexadecene,
            786.24,
            357.2/MM*1.0E3,
            3.0E-3,
            0.1423,
            MM,
            Tc,
            Pc,
            ω,
            RK_params,
            PR_params,
            PRSV_params,
            TTps
        )

    end


end



function rho(P,T,fluid::Hexadecene,EOS::PRSV)
    ## Based in the PRSV EOS
    v=solve_EOS(EOS,fluid,P,T)

    return fluid.MM/v
    # Tr=T/fluid.Tc
    # Av=0.808
    # Bv=2.575

    # Dv=-0.07844500194634374

    # v_new=(v+Dv)*(1-exp(Av-Bv/Tr^3))
    # return fluid.MM/v

    # return fluid.MM/v_new

end


function rho(P,T,fluid::Hexadecene,EOS::RK)

    v=solve_EOS(EOS,fluid,P,T)

    return fluid.MM/v
end


function rho(P,T,fluid::Hexadecene,EOS::PR)

    v=solve_EOS(EOS,fluid,P,T)

    return fluid.MM/v


    # return 0.5
end


function rho_stdP_hexadecene(T,hex::Hexadecene)
    return 4.608E2+8.18429E-1*(hex.Tc-T)-1.6747E-4(hex.Tc-T)^2
end

function rho(P,T,fluid::Hexadecene,EOS::TT)

    rho_stdP=rho_stdP_hexadecene(T,fluid)

    (E,F0,F1,F2)=fluid.TammannTait_params
    num=P-F0-F1*T-F2*T*T
    den=101325.0-F0-F1*T-F2*T*T

    return rho_stdP*(1-E*log(num/den))
end






rho(P,T,fluid::Hexadecene)=rho(P,T,fluid,TT())




function visc(P,T,fluid::Hexadecene)
    return (10^(10^9.975 * T^(-4.137))-0.8)*1E-3
end

function Cp(P,T,fluid::Hexadecene)
    # return fluid.cₚ
    return 1.38547E3+2.042*(T)+1.8283E-3*(T)^2
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
