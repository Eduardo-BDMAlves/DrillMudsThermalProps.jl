export Bentonite

struct Bentonite <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    function Bentonite()
        new(
            :Bentonite,
            2600.311,
            1150.0,
            1.3
        )
    end
end




function Cp(P,T,solid::Bentonite)
    return 12.187T^0.736304
end



function therm_cond(P,T,solid::Bentonite)
    # return solid.k
    return 0.0007699999999999954T-0.14974549999999853
end
