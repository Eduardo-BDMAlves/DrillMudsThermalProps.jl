module DrillMudsThermalProps


import Zygote: gradient
import Interpolations: LinearInterpolation

# Write your package code here.

include("Materials.jl")
include("Correlations.jl")

include("Mixture.jl")
include("ConductivityMixtures.jl")



end