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





function therm_cond(P, T, solid::Barite)
    Tc = T - 273.15
    return -9.0E-5 * Tc + 0.0974
end





