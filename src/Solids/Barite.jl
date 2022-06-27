export Barite

struct Barite <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    function Barite()
        new(
            :Barite,
            4100.0,
            4600.0,
            0.09515
        )
    end
end



function rho(P,T,solid::Barite)
    return solid.ρ/((1+0.8E-4(T-288.705))*(1-1E-10*(P-101325.0)))
end


function therm_cond(P, T, solid::Barite)
    Tc = T - 273.15
    return -9.0E-5 * Tc + 0.0974
end





