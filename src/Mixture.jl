export Mud,
        DrillFluid


abstract type Mud end


struct DrillFluid <: Mud
    Liquids::Vector{Fluid}
    Solids::Vector{Solid}

    nLiq::UInt
    nSol::UInt

    fS::Vector{Float64}
    fL::Vector{Float64}

    wt::Vector{Float64}


    function DrillFluid()
        #TODO: completar esse trecho com a primeira parte dos modelos de mistura... implementar mais de uma versão, um com o uso de frações de massa, outro com volume. Por enquanto não usar 

    end


end






