export Tsederberg,
    Maxwell,
    MaxwellMixture,
    Rayleigh,
    RayleighMixture,
    Churchill,
    ChurchillSolidComp,
    MaxwellSolidComp,
    RayleighSolidComp,
    Series,
    Parallel



## Model name

abstract type CondMixture end

struct Tsederberg <: CondMixture end

struct Maxwell <:CondMixture end

struct MaxwellMixture <: CondMixture end

struct MaxwellSolidComp <: CondMixture end

struct Rayleigh <: CondMixture
    a::Float64
    function Rayleigh(config::Symbol = :BCC)
        if config == :SCA
            new(
                1.31
            )
        elseif config == :BCC
            new(
                0.129
            )
        elseif config == :FCC
            new(
                0.0752
            )
        else
            new(
                0.129
            )
        end
    end
end

struct RayleighMixture <: CondMixture
    a::Float64
    function RayleighMixture(config::Symbol= :BCC)
        if config == :SCA
            new(
                1.31
            )
        elseif config == :BCC
            new(
                0.129
            )
        elseif config == :FCC
            new(
                0.0752
            )
        else
            new(
                0.129
            )
        end
    end
end

struct RayleighSolidComp <: CondMixture
    a::Float64
    function RayleighSolidComp(config::Symbol= :BCC)
        if config == :SCA
            new(
                1.31
            )
        elseif config == :BCC
            new(
                0.129
            )
        elseif config == :FCC
            new(
                0.0752
            )
        else
            new(
                0.129
            )
        end
    end
end

struct Churchill <: CondMixture end

struct ChurchillSolidComp <: CondMixture end

struct Series <: CondMixture end

struct Parallel <: CondMixture end

struct ChurchillLiqSolComp <: CondMixture end


## Functions

function therm_cond(P,T,Mud::DrillFluid, model::Tsederberg)
    return dot(therm_cond.(P,T,Mud.Liquids),Mud.wtL)+dot(therm_cond.(P,T,Mud.Solids),Mud.wtS)
end

function therm_cond(P,T,Mud::DrillFluid, model::Maxwell)

    fluid_list = Mud.Liquids
    solid_list = Mud.Solids

    (fl,fs) = update_fs(P,T,Mud)

    max_index=findfirst(fl.==maximum(fl))
    disc_index=findfirst(fs.==maximum(fs))

    k_c=therm_cond(P,T,fluid_list[max_index])

    k_d=therm_cond(P,T,solid_list[disc_index])
    f_d=fs[max_index]/(fl[max_index]+fs[disc_index])


    return k_c*(2.0k_c+k_d+2.0f_d*(k_d-k_c))/(2.0k_c+k_d-f_d*(k_d-k_c))
end

function therm_cond(P,T,Mud::DrillFluid, model::MaxwellMixture)
    
    (fl,fs) = update_fs(P,T,Mud)

    sum_fl=sum(fl)
    sum_fs=sum(fs)

    k_eq_L = dot(therm_cond.(P, T, Mud.Liquids), fl) / sum_fl
    k_eq_S = dot(therm_cond.(P, T, Mud.Solids), fs) / sum_fs

    return k_eq_L*(2.0k_eq_L+k_eq_S+2.0sum_fs*(k_eq_S-k_eq_L))/(2.0k_eq_L+k_eq_S-sum_fs*(k_eq_S-k_eq_L))
end


function therm_cond(P,T,Mud::DrillFluid, model::Rayleigh)

    fluid_list = Mud.Liquids
    solid_list = Mud.Solids

    (fl,fs) = update_fs(P,T,Mud)

    max_index=findfirst(fl.==maximum(fl))
    disc_index=findfirst(fs.==maximum(fs))

    k_c=therm_cond(P,T,fluid_list[max_index])

    k_d=therm_cond(P,T,solid_list[disc_index])
    f_d=fs[max_index]/(fl[max_index]+fs[disc_index])

    num_1=2.0+k_d/k_c
    num_2=1.0-k_d/k_c
    den_1=num_2
    den_2=4/3+k_d/k_c
    A=num_1/den_1
    B=num_2/den_2
    return k_c*(1.0-3.0f_d/(A+f_d-B*model.a*f_d^(10/3)))
end



function therm_cond(P,T,Mud::DrillFluid, model::RayleighMixture)

    (fl,fs) = update_fs(P,T,Mud)

    sum_fl=sum(fl)
    sum_fs=sum(fs)

    k_eq_L = dot(therm_cond.(P, T, Mud.Liquids), fl) / sum_fl
    k_eq_S = dot(therm_cond.(P, T, Mud.Solids), fs) / sum_fs

    num_1=2.0+k_eq_S/k_eq_L
    num_2=1.0-k_eq_S/k_eq_L
    den_1=num_2
    den_2=4/3+k_eq_S/k_eq_L
    A=num_1/den_1
    B=num_2/den_2
    return k_eq_L*(1.0-3.0sum_fs/(A+sum_fs-B*model.a*sum_fs^(10/3)))
end


function therm_cond(P,T,Mud::DrillFluid, model::Churchill)
    k_c=dot(therm_cond.(P,T,Mud.Liquids),Mud.wtL)

    (fl,fs) = update_fs(P,T,Mud)

    fS=sum(fs)

    k_d = dot(therm_cond.(P, T, Mud.Solids), fs) / fS

    lambda=k_d/k_c

    num_1=2.0+lambda
    num_2=6.0+3.0lambda
    num_3=3.0-3.0lambda

    den_1=1.0-lambda
    den_2=4.0+3.0lambda
    den_3=den_2

    fs73=fS^(7/3)
    fs103=fS^(10/3)

    return k_c*(num_1/den_1-2.0fS+0.409*num_2*fs73/den_2-2.133*num_3*fs103/den_3)/(num_1/den_1+fS+0.409num_2*fs73/den_2-0.906num_3*fs103/den_3)
end

function therm_cond(P,T,Mud::DrillFluid,model::ChurchillSolidComp)
    
    correct_wtL=Mud.wtL./sum(Mud.wtL)
    k_c=dot(therm_cond.(P,T,Mud.Liquids),correct_wtL)

    (fl,fs) = update_fs(P,T,Mud)

    fC=sum(fl)

    k_eff=k_c
    for (i,s) in enumerate(Mud.Solids)
        k_d=therm_cond(P,T,s)
        # fS=fs[i]
        fS=fs[i]/(fC+fs[i])
        k_c=k_eff

        lambda=k_d/k_c

        num_1=2.0+lambda
        num_2=6.0+3.0lambda
        num_3=3.0-3.0lambda

        den_1=1.0-lambda
        den_2=4.0+3.0lambda
        den_3=den_2

        fs73=fS^(7/3)
        fs103=fS^(10/3)

        k_eff = k_c*(num_1/den_1-2.0fS+0.409*num_2*fs73/den_2-2.133*num_3*fs103/den_3)/(num_1/den_1+fS+0.409num_2*fs73/den_2-0.906num_3*fs103/den_3)
        fC+=fs[i]
    end
    
    return k_eff
end



function therm_cond(P,T,Mud::DrillFluid, model::RayleighSolidComp)
    correct_wtL=Mud.wtL./sum(Mud.wtL)
    (fl,fs) = update_fs(P,T,Mud)

    fC=sum(fl)
    
    k_c=dot(therm_cond.(P,T,Mud.Liquids),correct_wtL)

    k_eff=k_c
    for (i,s) in enumerate(Mud.Solids)
        k_d=therm_cond(P,T,s)

        fS=fs[i]/(fC+fs[i])
        k_c=k_eff

        num_1=2.0+k_d/k_c
        num_2=1.0-k_d/k_c
        den_1=num_2
        den_2=4/3+k_d/k_c
        A=num_1/den_1
        B=num_2/den_2
        k_eff=k_c*(1.0-3.0fS/(A+fS-B*model.a*fS^(10/3)))

        fC+=fs[i]
    end
    return k_eff
end


function therm_cond(P,T,Mud::DrillFluid, model::MaxwellSolidComp)
    correct_wtL=Mud.wtL./sum(Mud.wtL)
    (fl,fs) = update_fs(P,T,Mud)

    fC=sum(fl)
    
    k_c=dot(therm_cond.(P,T,Mud.Liquids),correct_wtL)

    k_eff=k_c
    for (i,s) in enumerate(Mud.Solids)
        k_d=therm_cond(P,T,s)

        fS=fs[i]/(fC+fs[i])

        k_c=k_eff


        k_eff=k_c*(2.0k_c+k_d+2.0fS*(k_d-k_c))/(2.0k_c+k_d-fS*(k_d-k_c))
        fC+=fs[i]
    end
    return k_eff
end

function therm_cond(P,T,Mud::DrillFluid, model::Series)
    correct_wtL=Mud.wtL./sum(Mud.wtL)
    (fl,fs) = update_fs(P,T,Mud)

    fC=sum(fl)
    
    k_c=dot(therm_cond.(P,T,Mud.Liquids),correct_wtL)

    k_eff=k_c
    for (i,s) in enumerate(Mud.Solids)
        k_d=therm_cond(P,T,s)

        fS=fs[i]/(fC+fs[i])

        k_c=k_eff

        lamb=k_d/k_c

        k_eff=k_c*(lamb/(fS+(1-fS)*lamb))
        fC+=fs[i]
    end
    return k_eff
end

function therm_cond(P,T,Mud::DrillFluid, model::Parallel)
    correct_wtL=Mud.wtL./sum(Mud.wtL)
    (fl,fs) = update_fs(P,T,Mud)

    fC=sum(fl)
    
    k_c=dot(therm_cond.(P,T,Mud.Liquids),correct_wtL)

    k_eff=k_c
    for (i,s) in enumerate(Mud.Solids)
        k_d=therm_cond(P,T,s)

        fS=fs[i]/(fC+fs[i])

        k_c=k_eff

        lamb=k_d/k_c

        k_eff=k_c*((1-fS)+fS*lamb)
        fC+=fs[i]
    end
    return k_eff

end

function therm_cond(P,T,Mud::DrillFluid, model::ChurchillLiqSolComp)
    

    Full_list=list_sort_f(P,T,Mud)

    k_c=therm_cond(P,T,Full_list[1][2])
    fC=Full_list[1][1]


    k_eff=k_c
    for (fs,F) in Full_list[2:end]
        k_d=therm_cond(P,T,F)
        # fS=fs[i]
        fS=fs/(fC+fs)
        k_c=k_eff

        lambda=k_d/k_c

        num_1=2.0+lambda
        num_2=6.0+3.0lambda
        num_3=3.0-3.0lambda

        den_1=1.0-lambda
        den_2=4.0+3.0lambda
        den_3=den_2

        fs73=fS^(7/3)
        fs103=fS^(10/3)

        k_eff = k_c*(num_1/den_1-2.0fS+0.409*num_2*fs73/den_2-2.133*num_3*fs103/den_3)/(num_1/den_1+fS+0.409num_2*fs73/den_2-0.906num_3*fs103/den_3)
        fC+=fs
    end
    
    return k_eff
end





## General funcion

function therm_cond(P,T,Mud::DrillFluid) 
    return therm_cond(P,T,Mud,Tsederberg())
end







