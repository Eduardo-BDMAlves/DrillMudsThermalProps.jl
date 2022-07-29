export Barite

struct Barite <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    MM::Float64

    function Barite()

        MM=233.39
        Cp_mol=101.8
        new(
            :Barite,
            4100.0, #kg/m3
            Cp_mol/MM*1000, # cp
            0.09515, #k
            MM, #g/mol
        )
    end
end



function rho(P,T,solid::Barite)
    return solid.ρ/((1+0.8E-4(T-288.705))*(1-1E-10*(P-101325.0)))
end


function therm_cond(P, T, solid::Barite)
    # Tc = T - 273.15
    # return -9.0E-5 * Tc + 0.0974
    return 0.1219835-9E-5*T
end
