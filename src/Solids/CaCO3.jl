export CaCO3

struct CaCO3 <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    k::Float64

    # http://www.matweb.com/search/datasheet_print.aspx?matguid=bea4bfa9c8bd462093d50da5eebe78ac
    # https://thermtest.com/thermal-resources/materials-database
    # Caenn -> drilling mud
    function CaCO3()
        new(
            :CaCO3,
            2610.0,
            834.0,
            2.259
        )
    end
end