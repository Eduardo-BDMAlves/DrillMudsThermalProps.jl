export TestSolid,
    TestFluid,
    Fluid,
    Solid,
    Fluid_sets


abstract type Material end

abstract type Solid <: Material end
abstract type Fluid <: Material end

struct TestSolid <: Solid
    name::Symbol
    ρ::Float64
    cₚ::Float64
    μ::Float64
    k::Float64

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




## General computation common to every fluid

function (F::Fluid)(P, T)
    ρ = rho(P, T, F)
    κ = compress(P, T, F)
    β = therm_expand(P, T, F)
    C = Cp(P, T, F)
    k = therm_cond(P, T, F)
    μ = visc(P, T, F)

    return (ρ = ρ,
        κ = κ,
        β = β,
        Cₚ = C,
        k = k,
        μ = μ)
end
