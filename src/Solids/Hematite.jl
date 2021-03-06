export Hematite

struct Hematite <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    function Hematite()
        new(
            :Hematite,
            4900.0,
            649.8627065488718,
            6.413265500000001
        )
    end
end








function therm_cond(P, T, solid::Hematite)
    return (0.0839-6.63E-5*T)*100.0
end
