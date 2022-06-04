export nDodecane

struct nDodecane <: Fluid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    MM::Float64

    Tc::Float64
    Pc::Float64
    ω::Float64


    PR_params::NamedTuple

    function nDodecane()
        MM=170.33484
        # Acording to Lemmon2004
        Tc=658.1
        Pc=1.817E6
        # Acording to coolprop - http://www.coolprop.org/fluid_properties/fluids/n-Dodecane.html
        ω=0.574182221240689


        new(
            :nDodecane,
            750.0,
            373.30/MM*1000.0,
            0.00134,
            9999.0,
            MM,
            Tc,
            Pc,
            ω,
            set_up_EOS(PR(),Tc,Pc,ω)
        )



    end
end

rho(P,T,fluid::nDodecane) = rho(P,T,fluid,PR())


function rho(P,T,fluid::nDodecane,EOS::PR)

    v=solve_EOS(EOS,fluid,P,T)
    return fluid.MM/v
end








