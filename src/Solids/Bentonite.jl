export Bentonite

struct Bentonite <: Solid
    name::Symbol
    ρ::Float64
    Cₚ::Float64
    k::Float64

    function Bentonite()
        new(
            :Bentonite,
            2600.311,
            1150.0,#
            1.3#https://www.andra.fr/mini-sites/lille2007/abstract_lille2007/donnees/pdf/579_580_P_THME_19.pdf
        )
    end
end