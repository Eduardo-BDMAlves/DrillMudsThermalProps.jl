export Mud,
        DrillFluid,
        CTEVolDist,
        VarVolDist,
        update_fs,
        rho,
        Cp


abstract type Mud end


struct DrillFluid <: Mud
    Liquids::Vector{<:Fluid}
    Solids::Vector{<:Solid}

    nLiq::Integer
    nSol::Integer

    fL::Vector{Float64}
    sum_fL::Float64
    fS::Vector{Float64}
    sum_fS::Float64

    # wt::Vector{Float64}
    wtL::Vector{Float64}
    wtS::Vector{Float64}

    std_rho::Float64

    function DrillFluid(Flist::Vector{<:Fluid},Slist::Vector{<:Solid};fs::AbstractVector,normalize=false)

        nL = Integer(length(Flist))
        nS = Integer(length(Slist))

        missing_pos = findall(ismissing,fs)
        

        if !isempty(missing_pos)
            if length(missing_pos)>1
                @error "Too many volumetric proportions field left open. Package can complete the mixture with only one component at a time."
            end
            total_frac=sum(skipmissing(fs))
            if total_frac>1.0 
                @error "Normalization only available if no missing values are provided." 
            end

            fs[missing_pos[1]]=1.0-total_frac
            total_frac=1.0
        else
            total_frac=sum(fs)
            if total_frac != 1.0
                if normalize
                    @warn "Total volumetric fractions greater than 100%, normalizing volume to crrect. If that is not desired pass `normalize` as false."

                    fs=fs./total_frac
                    total_frac=1.0

                else
                    @error "Total volume fraction provided greater than 100%, if the values provided are total volumes change the keyword argument `normalize` to true and the values will be normalized."
                end
            end
        end

        fL=fs[1:length(Flist)]
        sum_fL=sum(fL)
        fS=fs[length(Flist)+1:end]
        sum_fS=sum(fS)

        Lrhos=[rho.ρ for rho ∈ Flist]
        Srhos=[rho.ρ for rho ∈ Slist]

        mLs=[rho*f for (rho,f) ∈ zip(Lrhos,fL)]
        mSs=[rho*f for (rho,f) ∈ zip(Srhos,fS)]

        mT=sum(mLs)+sum(mSs)

        wtL=mLs./mT
        wtS=mSs./mT

        std_rho=dot(fL,[rho.ρ for rho ∈ Flist])+dot(fS,[rho.ρ for rho ∈ Slist])

        new(
            Flist,
            Slist,
            nL,
            nS,
            fL,
            sum_fL,
            fS,
            sum_fS,
            wtL,
            wtS,
            std_rho
        )


    end


    # function DrillFluid(Flist,Slist;fs=missing,wts::Vector{Float64})
    #     #TODO: completar esse trecho com a primeira parte dos modelos de mistura... implementar mais de uma versão, um com o uso de frações de massa, outro com volume. Por enquanto não usar 
        



    # end



end







## Functions




function update_fs(P,T,fluid::DrillFluid)
    rho_T = rho(P,T,fluid)
    fl = fluid.wtL ./ rho.(P,T,fluid.Liquids).*rho_T

    fs = fluid.wtS ./ rho.(P,T,fluid.Solids).*rho_T

    return (fl,fs)
end





#for density

abstract type RhoModel end

struct CTEVolDist<:RhoModel end

struct VarVolDist<:RhoModel end

# density of discret phases

function rho(P,T,fluid::DrillFluid)
    return rho(P,T,fluid,VarVolDist())
end


function rho(P,T,fluid::DrillFluid,model::CTEVolDist)
    fluid_list=fluid.Liquids
    solids_list=fluid.Solids

    fL=fluid.fL
    fS=fluid.fS

    rho_l=[rho(P,T,f) for f ∈ fluid_list]
    rho_s=[rho(P,T,s) for s ∈ solids_list]

    rho_eq_l=sum(fL'*rho_l)
    rho_eq_s=sum(fS'*rho_s)

    f_TL=sum(fL)


    rho_T=rho_eq_l+rho_eq_s

    return rho_T
end


function rho(P,T,fluid::DrillFluid,model::VarVolDist)
    fluid_list=fluid.Liquids

    fL=fluid.fL

    rho_l=rho.(P,T,fluid.Liquids)

    rho_l_std=[f.ρ for f ∈ fluid_list]

    d_rho_l=rho_l-rho_l_std

    d_rho_l_dim=d_rho_l./rho_l

    return fluid.std_rho/(1.0-dot(fL,d_rho_l_dim))
end



abstract type CpModel end

struct CTEVolCp<:CpModel end

struct VarVolCp<:CpModel end
struct MassCp<:CpModel end

# density of discret phases

function Cp(P,T,fluid::DrillFluid)
    return Cp(P,T,fluid,MassCp())
end


function Cp(P,T,fluid::DrillFluid,mod::CTEVolCp)

    CpL=Cp.(P,T,fluid.Liquids)
    CpS=Cp.(P,T,fluid.Solids)

    return dot(fluid.fL,CpL)+dot(fluid.fS,CpS)
end

function Cp(P,T,fluid::DrillFluid,mod::VarVolCp)

    (fL,fS)=update_fs(P,T,fluid)

    CpL=Cp.(P,T,fluid.Liquids)
    CpS=Cp.(P,T,fluid.Solids)

    return dot(fL,CpL)+dot(fS,CpS)
end

function Cp(P,T,fluid::DrillFluid,mod::MassCp)

    CpsL=Cp.(P,T,fluid.Liquids)
    CpsS=Cp.(P,T,fluid.Solids)

    return dot(CpsL,fluid.wtL)+dot(CpsS,fluid.wtS)
end


