using DrillMudsThermalProps
using Documenter

DocMeta.setdocmeta!(DrillMudsThermalProps, :DocTestSetup, :(using DrillMudsThermalProps); recursive=true)

makedocs(;
    modules=[DrillMudsThermalProps],
    authors="Eduardo Alves",
    repo="https://github.com/Eduardo-BDMAlves/DrillMudsThermalProps.jl/blob/{commit}{path}#{line}",
    sitename="DrillMudsThermalProps.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Eduardo-BDMAlves.github.io/DrillMudsThermalProps.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Eduardo-BDMAlves/DrillMudsThermalProps.jl",
    devbranch="Development",
)
