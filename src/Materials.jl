
export TestSolid,
    TestFluid,
    Water,
    Barite


abstract type Material end

abstract type Solid <: Material end
abstract type Fluid <: Material end

# struct Water <: Fluid end

struct TestSolid{F} <: Solid
    ρ::F
    cₚ::F
    μ::F
    k::F
end


struct TestFluid <: Fluid
    name
    ρ
    cₚ
    μ
    k

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
            10000.0,
            0.949
        )
    end
end









