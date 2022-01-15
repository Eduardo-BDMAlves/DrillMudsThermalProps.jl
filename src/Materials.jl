
export TestSolid,
    TestFluid,
    Water,
    Water_std




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


struct TestFluids{F} <: Fluid
    ρ::F
    cₚ::F
    μ::F
    k::F
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














