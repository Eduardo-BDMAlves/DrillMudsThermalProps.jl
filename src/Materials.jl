
export TestSolid,
    TestFluid,
    Water


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

    return (ρ=ρ,missing)
end









## Solid

struct Barite <: Solid

    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

    function Water()
        new(
            :Barite,
            2000.0,
            5000.0,
            1.002E-3,
            0.568
        )
    end
end









