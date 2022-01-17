export rho,Cp,visc


## rho

function rho(P,T,fluid::Fluid)
    return missing
end

function rho(P,T,fluid::Water)
    T2=T*T
    PT=P*T
    PT2=P*T2
    a=8.79278807e+02
    b=1.11820207e+00
    c=-2.40938185e-03
    d=1.29027705e-06
    e=-5.46746114e-09
    f=8.06188264e-12

    return a+b*T+c*T2+d*P+e*PT+f*PT2
end



## Cp

function Cp(P,T,fluid::Fluid)
    return missing
end


function Cp(P,T,fluid::Water)
    b=3.70641927e+03
    c=7.13116317e-01
    d=1.05153132e-03
    e=8.59719978e-07
    f=-4.37409377e-09
    g=-4.26280988e-15
    h=1.14112829e-17

    P2=P*P
    
    PT=P*T

    return b + 2*c*T + e*P + 2*f*P*T + g*P2 + 2*h*P2*T
end




## k

function therm_cond(P,T,fluid::Fluid)
    return missing
end


function therm_cond(P,T,fluid::Water)
    a=-3.98200943e-01
    b= 5.32706765e-03
    c=-6.55323715e-06
    d= 1.19857114e-09
    e=-5.13478793e-12
    f= 9.20668097e-15
    g=-2.77585855e-24

    T2=T*T
    PT=P*T
    P2=P*P

    return a + b*T + c*T2 + d*P + e*PT + f*P*T2 + g*P2*T2
end



# mu


function visc(P,T,fluid::Fluid)
    return missing
end


function visc(P,T,fluid::Water)
    a = 2.65387609e+00
    b = 2.37561266e+02
    c =-2.37488893e+02
    d =-6.43059846e+08
    e = 6.43115832e+08
    f = 7.34303718e-07
    g =-9.76896404e-01

    return a*(exp((b-T)/(c+T)*((P+d)/(P-e)))) + f*T + g
end












