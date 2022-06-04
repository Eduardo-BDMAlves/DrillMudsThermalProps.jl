export Einstein


abstract type ViscMixture end

struct Einstein <: ViscMixture end
struct GuthSimha <: ViscMixture end

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

function visc(P,T,Mud::DrillFluid)
    return visc(P,T,Mud,GuthSimha())
end

