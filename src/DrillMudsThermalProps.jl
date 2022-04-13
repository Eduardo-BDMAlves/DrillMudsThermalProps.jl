module DrillMudsThermalProps


import Zygote: gradient
import Interpolations: LinearInterpolation
using JLD2
using Unitful
using DataFrames


# Write your package code here.

include("Erros.jl")
include("Materials.jl")
include("Correlations.jl")

include("Mixture.jl")
include("ConductivityMixtures.jl")




include("InputTableParser.jl")
include("Exporter.jl")



end