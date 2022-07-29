export Einstein


abstract type ViscMixture end

struct Einstein <: ViscMixture end
struct GuthSimha <: ViscMixture end
struct SingFluidGuthSimha <: ViscMixture end
struct SingFluidEinstein <: ViscMixture end

function visc(P,T,Mud::DrillFluid, model::Einstein)
    # cont_index=findfirst(Mud.fL .== maximum(Mud.fL))
    (fL,fS)=update_fs(P,T,Mud)
    cont_index=findfirst(fL .== maximum(fL))
    return visc(P,T,Mud.Liquids[cont_index])*(1.0+2.5*sum(fS))
end

function visc(P,T,Mud::DrillFluid, model::GuthSimha)
    (fL,fS)=update_fs(P,T,Mud)
    cont_index=findfirst(fL .== maximum(fL))
    return visc(P,T,Mud.Liquids[cont_index])*(1.0+2.5*sum(fS)+14.1*sum(fS)^2)
end


function visc(P,T,Mud::DrillFluid, model::SingFluidGuthSimha)
    (fL,fS)=update_fs(P,T,Mud)
    tfL=sum(fL)
    tfS=sum(fS)

    if length(fL) > 2
        return NaN
    elseif length(fL)==1
        mu=visc(P,T,Mud.Liquids[1])
    else
        nfL= fL./tfL
        muexp=[visc(P,T,L)^(nfL[i]) for (i,L) in enumerate(Mud.Liquids)]
        mu=prod(muexp)
    end

    return mu*(1.0+2.5tfS+14.1tfS^2)
end

function visc(P,T,Mud::DrillFluid, model::SingFluidEinstein)
    (fL,fS)=update_fs(P,T,Mud)
    tfL=sum(fL)
    tfS=sum(fS)

    if length(fL) > 2
        return NaN
    elseif length(fL)==1
        mu=visc(P,T,Mud.Liquids[1])
    else
        nfL= fL./tfL
        muexp=[visc(P,T,L)^(nfL[i]) for (i,L) in enumerate(Mud.Liquids)]
        mu=prod(muexp)
    end

    return mu*(1.0+2.5tfS)
end


function visc(P,T,Mud::DrillFluid)
    return visc(P,T,Mud,SingFluidGuthSimha())
end

