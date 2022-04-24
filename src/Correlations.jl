export rho,
    Cp,
    visc,
    therm_cond,
    therm_expand,
    compress,
    std_props

## Std props

function std_props(fluid::Fluid)
    ρ = fluid.ρ
end


## rho

function rho(P, T, fluid::Fluid)
    @error "Property not yet implemented for the fluid. - ρ"
    throw(NotImplementedError())
    return missing
end


## Cp

function Cp(P, T, fluid::Fluid)
    @error "Property not yet implemented for the fluid. - Cₚ"
    throw(NotImplementedError())
    return missing
end


## k

function therm_cond(P, T, fluid::Fluid)
    @error "Property not yet implemented for the fluid. - k"
    throw(NotImplementedError())
    return missing
end


# mu


function visc(P, T, fluid::Fluid)
    @error "Property not yet implemented for the fluid. - μ"
    throw(NotImplementedError())
    return missing
end





# therm expand

function therm_expand(P, T, fluid::Fluid)
    return -gradient(x -> rho(P, x, fluid), T)[1] / rho(P, T, fluid)
end


## compressibility

function compress(P, T, fluid::Fluid)
    return gradient(x -> rho(x, T, fluid), P)[1] / rho(P, T, fluid)
end



# Solid

## density

function rho(P, T, solid::Solid)
    # @error "Property not yet implemented for solid - ρ"
    # throw(NotImplementedError())
    # return missing
    return solid.ρ
end

## Cp

function Cp(P, T, solid::Solid)
    # @error "Property not yet implemented for solid - Cₚ"
    # throw(NotImplementedError())
    # return missing
    return solid.Cₚ
end

function Cp(P, T, solid::Barite)
    return 4600.0
end

function Cp(P, T, solid::CaCO3)
    return solid.cₚ
end

function Cp(P,T,solid::Hematite)
    A=93.43834
    B=108.3577
    C=-50.86447
    D=25.58683
    E=-1.611330
    t=T/1000.0
    CpM=A+B*t+C*t^2+D*t^3+E/(t^2)
    return CpM/159.688*1000.0
end

## therm cond

function therm_cond(P, T, solid::Solid)
    @error "Property not yet implemented for solid - k"
    throw(NotImplementedError())
    return missing
end



# function therm_cond(P, T, solid::CaCO3)
#     return solid.k
# end




