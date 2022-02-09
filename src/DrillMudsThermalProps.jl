module DrillMudsThermalProps


using Zygote
using Interpolations

# Write your package code here.

include("Materials.jl")
include("Correlations.jl")

include("Mixture.jl")
include("ConductivityMixtures.jl")



end
