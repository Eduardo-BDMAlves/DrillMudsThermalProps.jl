
export TestSolid,
    TestFluid,
    Water,
    Barite,
    Fluid,
    Solid,
    Fluid_sets,
    Brine_NaCl


abstract type Material end

abstract type Solid <: Material end
abstract type Fluid <: Material end

# Fluid_sets=Union{(Fluid, (Fluid)...)...}

# struct Water <: Fluid end

struct TestSolid <: Solid
    name
    ρ
    cₚ
    μ
    k

    function TestSolid()
        new(
            :TestSolid,
            1000.0,
            4000.0,
            1.002E-3,
            0.58
        )
    end
end


struct TestFluid <: Fluid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    function TestFluid()
        new(
            :Water,
            1000.0,
            4000.0,
            1.002E-3,
            0.58
        )
    end
end




struct Water <: Fluid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    function Water()
        new(
            :Water,
            1000.0,
            4000.0,
            1.002E-3,
            0.58
        )
    end
end


struct Brine_NaCl <: Fluid

    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64
    x::Float64



    function BrineNaCl(ppg=9.8)

        ppgs=[8.6:0.1:9.8]
        xs=[
            0.45501,0.061330,0.077159,0.092986,0.108814,0.124642,0.140470,0.156297,0.172125,0.187953,0.203782,0.219610,0.235437
        ]

        itp=LinearInterpolation(ppgs,xs)

        x=itp(ppg)


        rho_std=0.0
        cp_std=0.0
        mu_std=0.0
        k_std=0.0

        new(
            :Brine_NaCl,
            rho_std,
            cp_std,
            mu_std,
            k_std,
            x
        )
    end


end










## General computation common to every fluid

function (F::Fluid)(P,T)
    ρ=rho(P,T,F)
    κ=compress(P,T,F)
    β=therm_expand(P,T,F)
    C=Cp(P,T,F)
    k=therm_cond(P,T,F)
    μ=visc(P,T,F)

    println("hello")
    return (ρ=ρ,
            κ=κ,
            β=β,
            Cₚ=C,
            k=k,
            μ=μ)
end









## Solid

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
            0.949
        )
    end
end









