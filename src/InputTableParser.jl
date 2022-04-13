export solicitation_table_parser



function solicitation_table_parser(df,target_rho,total_vol)

    lbbbl_to_kgm3=2.857142857142857u"kg/m^3"

    ppg_to_kgm3=119.83u"kg/m^3"

    galbbl_to_mLmL=0.0238

    mass=Float64[]
    mass_unit=Unitful.Quantity[]

    QSP_indexes=Integer[]
    vol_indexes=Integer[]

    water_index=Integer[]
    water_content=Fluid[]
    oil_index=Integer[]
    oil_content=Fluid[]
    solids_index=Integer[]
    solid_content=Solid[]

    aditive_indexes=Integer[]

    has_QSP=false
    has_pure_water=false
    pure_water_index=0
    len_components=size(df)[1]


    #check if QSP of water and salt
    comp_df=copy(df)
    if any(df.Componente .== :H2O) && any(df.Componente .== :NaCl)
        w_i=findall(df.Componente .== :H2O)
        Na_i=findall(df.Componente .== :NaCl)
        
        w_QSP=any()
        
    end


    for i in range(1,len_components)

        row=df[i,:]

        # println(row.Componente)
        if row.unit == "lb/bbl"
            convert=lbbbl_to_kgm3
            m=upreferred(row.conc*convert*total_vol)
            if row.Componente == :Obiturante
                append!(solid_content,[CaCO3()])
                append!(solids_index,i)
            elseif row.Componente == :Reologico
                ben=Bentonite()
                append!(solid_content,[ben])
                append!(solids_index,i)
            end
        elseif row.unit == "bbl/bbl" || row.unit == "v/v"
            vol=row.conc*total_vol
            append!(vol_indexes,i)
            rho=1.0
            if row.Componente == :Olefina || row.Componente == :Hexadeceno
                fluid=Hexadecene()
                rho=0.728*1000.0
                # rho=fluid.ρ
                append!(oil_index,i)
                append!(oil_content,[fluid])
            elseif row.Componente == :Brine_NaCl && !ismissing(row.conc)
                fluid=BrineNaCl(9.8)
                rho=fluid.ρ
                append!(water_index,i)
                append!(water_content,[fluid])
            elseif row.Componente == :Obiturante
                append!(solid_content,[CaCO3()])
                rho=CaCO3().ρ
                append!(solids_index,i)
            end

            m=upreferred(vol*rho*u"kg/m^3")
        elseif row.unit == "% v/v"
            vol=row.conc*total_vol/100

            if row.Componente == :Liovac4260
                rho=700.0
            end
            println("Liovac")
            m=upreferred(vol*rho*u"kg/m^3")
        elseif row.unit == "gal/bbl"
            vol=galbbl_to_mLmL*total_vol

            if row.Componente == :Polifoan
                rho=1.01u"g/cm^3"
            end

            m=upreferred(vol*rho)
        elseif ismissing(row.conc) && occursin(row.unit,"QSP")
            has_QSP=true
            append!(QSP_indexes,i)
            m=0.0u"kg"

            if row.Componente == :Barita
                append!(solids_index,i)
                append!(solid_content,[Barite()])
            elseif row.Componente == :H2O
                has_pure_water=true
                pure_water_index=i
            end




        end
        println((i=i,))
        append!(mass_unit,m)
        append!(mass,ustrip(m))
    end
    println(mass_unit)

    ## Encontra as massas dos QSP...
    fluid_rho=target_rho*ppg_to_kgm3
    total_mass=upreferred(fluid_rho*total_vol)
    if !isempty(QSP_indexes)
        if length(QSP_indexes) == 1
            m_a=0.0u"kg"
            m_B=0.0u"kg"
            for i in range(1,len_components)
                if df[i,:].aditive
                    append!(aditive_indexes,i)
                    m_a+=mass_unit[i]
                end
                if i in vol_indexes
                    m_B+=mass_unit[i]
                end
            end
            m_s=0.0u"kg"
            for j in solids_index
                if j ∉ QSP_indexes
                    m_s+=mass_unit[j]
                end
            end
            mass_unit[QSP_indexes[1]]=total_mass-m_a-m_B
        else 
            # if has_pure_water && any(findall(.Componente .== :NaCl))
            #     println("hello")
            #     @error "Not Implemented"
            # end
        end
    end

    ## Encontrar as frações volumétricas considerando aditivos dissolvidos na água.

    fP=[0.0,0.0,0.0]
    fcomp=[
        [],
        [],
        []
    ]

    vol_sol=Unitful.Quantity[]
    vol_oil=0.0u"m^3"
    vol_water=0.0u"m^3"


    for (i,f) in zip(oil_index,oil_content)
        rho=f.ρ*u"kg/m^3"
        vol=upreferred(mass_unit[i]/rho)
        vol_oil+=vol
        append!(fcomp[2],vol_oil)
    end

    for (i,s) in zip(solids_index,solid_content)
        rho=s.ρ*u"kg/m^3"
        vol=upreferred(mass_unit[i]/rho)
        append!(vol_sol,vol)
        append!(fcomp[3],vol)
    end




    vol_water=total_vol-sum(vol_oil)-sum(vol_sol)
    
    fP[1]=upreferred(vol_water/total_vol)
    fP[2]=upreferred(sum(vol_oil)/total_vol)
    fP[3]=upreferred(sum(vol_sol)/total_vol)

    fcomp[1] = [upreferred(vol_water/total_vol)]
    fcomp[2] .= upreferred.(fcomp[2]./total_vol)
    fcomp[3] = upreferred.(fcomp[3]./total_vol)
    
    mass_frac=mass_unit./total_mass



    
    df.mass=mass_unit
    (
        total_mass=total_mass,
        total_mass_no_unit=ustrip(total_mass),
        total_vol=total_vol,
        total_vol_no_unit=ustrip(upreferred(total_vol)),
        components_mass=mass_unit,
        components_mass_no_unit=mass,
        aditive_mass=m_a,
        mass_frac=mass_frac,
        fP=fP,
        fcomp=fcomp,
        rho_std=target_rho*ppg_to_kgm3,
        solid_indexes=solids_index,
        aditive_indexes=aditive_indexes,
        df=df,
        # sol_indexes=solids_index
    )



end









