export Tsederberg



## Model name

abstract type CondMixture end

struct Tsederberg <: CondMixture end

struct Maxwell <:CondMixture end

struct Rayleigh <: CondMixture
    alpha::Float64
end

struct Churcill <: CondMixture end


struct MaxwellMixtures <: CondMixture end










## Functions

function therm_cond(P,T,Mud::DrillFluid, model::Tsederberg)
    return dot(therm_cond.(P,T,Mud.Liquids),Mud.wtL)+dot(therm_cond.(P,T,Mud.Solids),Mud.wtS)
end


function therm_cond(P,T,Mud::DrillFluid, model::MaxwellMixtures)
    
end


function therm_cond(P,T,Mud::DrillFluid) 
    return therm_cond(P,T,Mud,Tsederberg())
end







