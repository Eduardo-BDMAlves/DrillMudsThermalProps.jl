export EOS,
PRSV,
PR,
RK,
TT


abstract type EOS end

abstract type EOS_cubic <:EOS end


struct PRSV <:EOS_cubic end
struct PR <:EOS_cubic  end
struct RK <: EOS_cubic end
struct TT <:EOS end


function set_up_EOS(EOS::PR,Tc,Pc,ω)
        
    R=8.31446261815324*1.0E3

    ac=0.45724*(R*Tc)^2/Pc

    b=0.07780*R*Tc/Pc

    return (
        ac=ac,
        b=b,
        ω=ω
    )
end

function set_up_EOS(EOS::PRSV,Tc,Pc,k0,k1)
    R=8.31446261815324*1.0E3

    ac=0.45724*(R*Tc)^2/Pc

    b=0.07780*R*Tc/Pc

    return (
            ac=ac,
            b=b,
            k0=k0,
            k1=k1
    )
end

function set_up_EOS(EOS::RK,Tc,Pc,ω)

    R=8.31446261815324*1.0E3

    ac=0.45724*(R*Tc)^2/Pc

    b=0.07780*R*Tc/Pc

    return (
        ac=ac,
        b=b,
        ω=ω
    )
end

function solve_cubic(C1,C2,C3)
    C4=C1/3

    D=0.5*(C2*C4-C3)-C4^3
    E=C4^2-C2/3
    Δ=D^2-E^3

    tol=1.0E-8
    v=0.0
    if Δ>tol

        F=cbrt(D+sqrt(Δ))
        G=cbrt(D-sqrt(Δ))
        v=F+G-C4

    elseif Δ≥-tol
        v1=2cbrt(D)-C4
        v2=-cbrt(D)-C4
        println((v1=v1,v2=v2))
        if v1 ≥ 0 && v2 ≥ 0
            v=min(v1,v2)
        elseif v1 ≥ 0
            v=v1 
        elseif v2 ≥ 0
            v=v2
        else
            v=Inf64
        end
    else
        J=acos(D/sqrt(E^3))
        v1=2sqrt(E)*cos(J/3)-C4
        v2=2sqrt(E)*cos(J/3 + 2π/3)-C4
        v3=2sqrt(E)*cos(J/3+4π/3)-C4


        if v1 ≥ 0 && v2 ≥ 0 && v3 ≥ 0
            v=min(v1,v2,v3)
        elseif v1 ≥ 0 && v2 ≥ 0
            v=min(v1,v2)
        elseif v2 ≥ 0 && v3 ≥ 0
            v=min(v2,v3)
        elseif v1 ≥ 0 && v3 ≥ 0
            v=min(v1,v3)
        elseif v1 ≥ 0
            v=v1 
        elseif v2 ≥ 0
            v=v2
        elseif v3 ≥ 0
            v=v3
        else
            v=Inf64
        end
    end


    return v
end


function solve_EOS(EOS::RK,fluid::Fluid,P,T)
    R=8.31446261815324E3

    a,b,ω=fluid.RK_params

    Tr=T/fluid.Tc

    alpha=(1+(0.37464+1.54226ω-0.26992ω^2)*(1-sqrt(Tr)))^2

    A=a*alpha*P/(R*T)^2
    B=b*P/(R*T)

    C1=(B-1)
    C2=A-2B-3B^2
    C3=(B^2+B-A)*B


    v=solve_cubic(C1,C2,C3)
    return v
end



function solve_EOS(EOS::PR,fluid::Fluid,P,T)

    R=8.31446261815324E3

    ac,b,ω=fluid.PR_params

    Tr=T/fluid.Tc

    alpha=(1+(0.37464+1.54226ω-0.26992*(ω^2))*(1-sqrt(Tr)))^2

    a=alpha*ac
#     println((α=alpha,ac=ac,b=b))

    C1=b-R*T/P
    C2=a/P-2b*R*T/P-3b*b
    # C3=b*(b*b+b-a/P_loc)
    C3=b*(b*b+b*R*T/P-a/P)

    v=solve_cubic(C1,C2,C3)
end


function solve_EOS(EOS::PRSV,fluid::Fluid,P,T)

    R=8.31446261815324*1.0E3

    ac,b,k0,k1=fluid.PRSV_params

    Tr=T/fluid.Tc
    κ=k0+k1*(1+sqrt(Tr))*(0.7-Tr)
    alpha=(1+κ*(1-sqrt(Tr)))*(1+κ*(1-sqrt(Tr)))

    C1=b-R*T/P
    C2=ac*alpha/P-2R*T*b/P-3b^2
    C3=(b^2+R*T*b/P-ac*alpha/P)*b

    v=solve_cubic(C1,C2,C3)
    return v

end