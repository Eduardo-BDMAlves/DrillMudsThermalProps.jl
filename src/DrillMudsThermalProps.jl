module DrillMudsThermalProps


import Zygote: gradient
import Interpolations: LinearInterpolation
using JLD2

# Write your package code here.

include("Erros.jl")
include("Materials.jl")
include("Correlations.jl")

include("Mixture.jl")
include("ConductivityMixtures.jl")





include("Exporter.jl")



end