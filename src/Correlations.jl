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

function rho(P, T, fluid::Water)
    a = 8.79278807E+02
    b = 1.11820207E+00
    c = -2.40938185E-03
    d = 1.29027705E-06
    e = -5.46746114E-09
    f = 8.06188264E-12

    return a + b * T + c * T^2 + d * P + e * P * T + f * P * T^2
end


function rho(P, T, fluid::BrineNaCl)

    x = fluid.x

    a, b, c, d, g, h, m, n, p, q, r, s = 8.23073018e+02, 1.38381011e+00, -2.85110020e-03, -9.11569454e-18, 4.04861229e-08, 2.38570427e-15, -2.65736884e-09, 9.52410475e-12, 1.27087646e+03, -7.83633987e-01, 6.26293103e-09, -2.11977575e-03

    d1 = a + b * T + c * T^2
    d2 = g * P + h * P^2
    d3 = m * P * T + n * P * T^2 + d * P^2 * T
    d4 = p * x + q * x * T + r * x * P * T + s * x * T^2

    return d1 + d2 + d3 + d4
end

## Cp

function Cp(P, T, fluid::Fluid)
    @error "Property not yet implemented for the fluid. - Cₚ"
    throw(NotImplementedError())
    return missing
end


function Cp(P, T, fluid::Water)
    b = 3.70641927E+03
    c = 7.13116317E-01
    d = 1.05153132E-03
    e = 8.59719978E-07
    f = -4.37409377E-09
    g = -4.26280988E-15
    h = 1.14112829E-17

    return b + 2 * c * T + e * P + 2 * f * P * T + g * P * 2 + 2 * h * P^2 * T
end

function Cp(P, T, fluid::BrineNaCl)

    x = fluid.x

    b, c, d, e, f, g, i, j, k, l, m = (2.48218846e+03, 3.80939153e-06, -3.45541918e-15, -2.57191794e+03, 3.85762099e+03, 1.93733545e-07, 2.27236796e+00, -9.07624069e-09, 1.27273598e-17, 1.38313250e+00, -1.96832796e+00)

    h_1 = b + c * P + d * P^2 + e * x + f * x^2 + g * x * P
    h_2 = 2 * (i + j * P + k * P^2 + l * x + m * x^2) * T
    return h_1 + h_2

end

## k

function therm_cond(P, T, fluid::Fluid)
    @error "Property not yet implemented for the fluid. - k"
    throw(NotImplementedError())
    return missing
end


function therm_cond(P, T, fluid::Water)
    a = -3.98200943E-01
    b = 5.32706765E-03
    c = -6.55323715E-06
    d = 1.19857114E-09
    e = -5.13478793E-12
    f = 9.20668097E-15
    g = -2.77585855E-24

    return a + b * T + c * T^2 + d * P + e * P * T + f * P * T^2 + g * P^2 * T^2
end


function therm_cond(P, T, fluid::BrineNaCl)

    x = fluid.x

    a, b, c, d, g, h, m, n, p, q, r, s, t = (-9.60221679e-01, 1.06170815e-02, -2.25969876e-05, 1.57190072e-08, -8.34381967e-10, -2.63526539e-19, 4.75581304e-12, -3.22756205e-15, 4.47716476e-01, -4.57135183e-03, 1.55078050e-09, -2.82355304e-12, 8.68547458e-06)

    k_1 = a + b * T + c * T^2 + d * T^3
    k_2 = g * P + h * P^2
    k_3 = m * P * T + n * P * T^2
    k_4 = p * x + q * x * T + r * x * P + s * x * P * T + t * x * T^2

    return k_1 + k_2 + k_3 + k_4
end


# mu


function visc(P, T, fluid::Fluid)
    @error "Property not yet implemented for the fluid. - μ"
    throw(NotImplementedError())
    return missing
end


function visc(P, T, fluid::Water)
    a = 2.65387609E+00
    b = 2.37561266E+02
    c = -2.37488893E+02
    d = -6.43059846E+08
    e = 6.43115832E+08
    f = 7.34303718E-07
    g = -9.76896404E-01

    return a * (exp((b - T) / (c + T) * ((P + d) / (P - e)))) + f * T + g
end


function visc(P, T, fluid::BrineNaCl)

    x = fluid.x

    a, b, c, d, g, h, m, p, q, r = (-9.29821669e+00, -5.62181704e+02, -2.60785572e-07, 2.24278114e+03, -1.75402843e-16, 3.70625444e+05, -4.33694111e+05, 4.47607787e-12, -1.44053809e-20, -1.33905655e-05)


    v_1 = a + (b + c * P + d * x + g * P^2) / T
    v_2 = (h + m * x) / T^2
    v_3 = p * P * T + q * P^2 * T + r * T^2 * x^2
    return exp(v_1 + v_2 + v_3)
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

function rho(P, T, fluid::Solid)
    @error "Property not yet implemented for solid - ρ"
    throw(NotImplementedError())
    return missing
end

function rho(P, T, solid::Barite)
    return solid.ρ
end

function rho(P, T, solid::CaCO3)
    return solid.ρ
end

function rho(P,T,solid::Hematite)
    return solid.ρ
end

## Cp

function Cp(P, T, solid::Solid)
    @error "Property not yet implemented for solid - Cₚ"
    throw(NotImplementedError())
    return missing
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

function therm_cond(P, T, solid::Barite)
    Tc = T - 273.15
    return -9.0E-5 * Tc + 0.0974
end


function therm_cond(P, T, solid::CaCO3)
    return solid.k
end


function therm_cond(P, T, solid::Hematite)
    return (0.0839-6.63E-5*T)*100.0
end


