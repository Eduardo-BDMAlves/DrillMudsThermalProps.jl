export rho,
        Cp,
        visc,
        therm_cond,
        therm_expand,
        compress,
        std_props

## Std props

function std_props(fluid::Fluid)
    ρ=fluid.ρ
end


## rho

function rho(P,T,fluid::Fluid)
    @error "Property not yet implemented for the fluid. - ρ"
    return missing
end

function rho(P,T,fluid::Water)
    a=8.79278807E+02
    b=1.11820207E+00
    c=-2.40938185E-03
    d=1.29027705E-06
    e=-5.46746114E-09
    f=8.06188264E-12

    return a+b*T+c*T^2+d*P+e*P*T+f*P*T^2
end



## Cp

function Cp(P,T,fluid::Fluid)
    @error "Property not yet implemented for the fluid. - Cₚ"
    return missing
end


function Cp(P,T,fluid::Water)
    b=3.70641927E+03
    c=7.13116317E-01
    d=1.05153132E-03
    e=8.59719978E-07
    f=-4.37409377E-09
    g=-4.26280988E-15
    h=1.14112829E-17

    return b + 2*c*T + e*P + 2*f*P*T + g*P*2 + 2*h*P^2*T
end


## k

function therm_cond(P,T,fluid::Fluid)
    @error "Property not yet implemented for the fluid. - k"
    return missing
end


function therm_cond(P,T,fluid::Water)
    a=-3.98200943E-01
    b= 5.32706765E-03
    c=-6.55323715E-06
    d= 1.19857114E-09
    e=-5.13478793E-12
    f= 9.20668097E-15
    g=-2.77585855E-24

    return a + b*T + c*T^2 + d*P + e*P*T + f*P*T^2 + g*P^2*T^2
end


# mu


function visc(P,T,fluid::Fluid)
    @error "Property not yet implemented for the fluid. - μ"
    return missing
end


function visc(P,T,fluid::Water)
    a = 2.65387609E+00
    b = 2.37561266E+02
    c =-2.37488893E+02
    d =-6.43059846E+08
    e = 6.43115832E+08
    f = 7.34303718E-07
    g =-9.76896404E-01

    return a*(exp((b-T)/(c+T)*((P+d)/(P-e)))) + f*T + g
end


# therm expand

function therm_expand(P,T,fluid::Fluid)
    return -gradient(x->rho(P,x,fluid),T)[1]/rho(P,T,fluid)
end


## compressibility

function compress(P,T,fluid::Fluid)
    return gradient(x->rho(x,T,fluid),P)[1]/rho(P,T,fluid)
end



# Solid

## density

function rho(P,T,fluid::Solid)
    @error "Property not yet implemented for solid - ρ"
    return missing
end

function rho(P,T,solid::Barite)
    return solid.ρ
end


## Cp

function Cp(P,T,solid::Solid)
    @error "Property not yet implemented for solid - Cₚ"
    return missing
end

function Cp(P,T,solid::Barite)
   return 4600.0
end


## therm cond

function therm_cond(P,T,solid::Solid)
    @error "Property not yet implemented for solid - k"
    return missing
end

function therm_cond(P,T,solid::Barite)
    Tc=T-273.15
    return -9.0E-5*Tc+0.0974
end






