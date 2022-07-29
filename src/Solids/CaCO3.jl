export CaCO3

struct CaCO3 <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64
    MM::Float64

    # http://www.matweb.com/search/datasheet_print.aspx?matguid=bea4bfa9c8bd462093d50da5eebe78ac
    # https://thermtest.com/thermal-resources/materials-database
    # Caenn -> drilling mud
    function CaCO3()
        MM=100.087
        Cp_mol=83.5


        new(
            :CaCO3,
            2710.0,
            Cp_mol/MM*1000,
            2.259,
            MM
        )
    end
end



function therm_cond(P,T,solid::CaCO3)
    return solid.k
end





