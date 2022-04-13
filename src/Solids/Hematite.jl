export Hematite

struct Hematite <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C1317608&Type=JANAFS&Table=on
    # https://www.researchgate.net/publication/252121179_Thermal_Conductivity_of_Magnetite_and_Hematite
    function Hematite()
        new(
            :Hematite,
            4900.0,
            649.8627065488718,
            6.413265500000001
        )
    end
end