export Tsederberg





abstract type CondMixture end

struct Tsederberg <: CondMixture end

function therm_cond(P,T,Mud::DrillFluid, model::Tsederberg)
    return dot(therm_cond.(P,T,Mud.Liquids),Mud.wtL)+dot(therm_cond.(P,T,Mud.Solids),Mud.wtS)
end


function therm_cond(P,T,Mud::DrillFluid) 
    return therm_cond(P,T,Mud,Tsederberg())
end







