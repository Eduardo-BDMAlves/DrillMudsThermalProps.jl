module DrillMudsThermalProps


import Zygote: gradient
import Interpolations: LinearInterpolation
using JLD2
using Unitful
using DataFrames


# Write your package code here.

include("Erros.jl")
include("Materials.jl")


Fluids_paths=joinpath(@__DIR__,"Fluids")
Fluids_files=joinpath.(Fluids_paths,readdir(Fluids_paths))
include.(Fluids_files)

Solids_paths=joinpath(@__DIR__,"Solids")
Solids_files=joinpath.(Solids_paths,readdir(Solids_paths))
include.(Solids_files)



include("Correlations.jl")

include("Mixture.jl")
include("ConductivityMixtures.jl")




include("InputTableParser.jl")
include("Exporter.jl")



end