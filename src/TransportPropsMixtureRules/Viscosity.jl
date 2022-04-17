export Einstein


abstract type ViscMixture end

struct Einstein <: ViscMixture end

function visc(P,T,Mud::DrillFluid, model::Einstein)
    cont_index=findfirst(Mud.fL .== maximum(Mud.fL))
    return visc(P,T,Mud.Liquids[cont_index])*(1.0+2.5*sum(Mud.fS))
end


function visc(P,T,Mud::DrillFluid) 
    return visc(P,T,Mud,Einstein())
end

