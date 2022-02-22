export generate_JLD2



function generate_JLD2(fname,props_dics,Ts,Ps,fluid;include_std=true)

    pair_PT_table=[(P=P,T=T) for P in Ps, T in Ts]
    pair_PT=[pair_PT_table...]

    len_PT=length(pair_PT)
    jump_T=length(Ts)
    size_PT=size(pair_PT_table)

    props_keys=keys(props_dics)



    props_lines=Dict{Union{Symbol,String},Array{Float64}}(
        [k => zeros(Float64,len_PT) for k in props_keys]
    )

    props_tabs=Dict{Union{Symbol,String},Array{Float64}}(
        [k => zeros(Float64,size_PT) for k in props_keys]
    )
 


    for k in props_keys
        for (i, P) in enumerate(Ps)
            for (j, T) in enumerate(Ts)
                prop_calc = props_dics[k](P, T, fluid)
    
                props_tabs[k][i, j] = prop_calc
    
                props_lines[k][(i-1)*jump_T+j] = prop_calc
    
            end
        end
    end




    jldopen(fname,"w",compress = true) do file
        file["fluid"]=fluid

        file["inputs/Ts"]=Ts
        file["inputs/Ps"]=Ps
        file["inputs/pairs"]=pair_PT
        file["inputs/table"]=pair_PT_table

        file["inputs/len_data"]=len_PT
        file["props"]=[k for k in props_keys]

        file["data/linear"]=props_lines
        file["data/tabular"]=props_tabs


    end


    return true

end




