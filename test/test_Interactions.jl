using Test, CSV
include("../src/basemodel.jl")
include("../src/components/Interactions.jl")

benchmark = CSV.read("../data/benchmark/Interactions.csv", DataFrame)

## Setup the model

model = test_model()
interact = addInteractions(model)

interact[:I_AMOC] = benchmark.I_AMOC
interact[:VGIS] = 1 .- benchmark."1 - V_GIS"
interact[:p_WAIS] = benchmark.I_WAIS
interact[:I_AMAZ] = benchmark.I_AMAZ
interact[:mNINO3pt4] = 1 .- benchmark."1 - m_NINO/m_0" ./ 10000

run(model)

## Test the model

f_AMOC = model[:Interactions,  :f_AMOC]
f_AMOC_compare = benchmark.f_AMOC
@test f_AMOC ≈ f_AMOC_compare rtol=1e-4

f_GIS = model[:Interactions,  :f_GIS]
f_GIS_compare = benchmark.f_GIS
@test f_GIS ≈ f_GIS_compare rtol=1e-4

f_WAIS = model[:Interactions,  :f_WAIS]
f_WAIS_compare = benchmark.f_WAIS
@test f_WAIS ≈ f_WAIS_compare rtol=1e-4

f_AMAZ = model[:Interactions,  :f_AMAZ]
f_AMAZ_compare = benchmark.f_AMAZ
@test f_AMAZ ≈ f_AMAZ_compare rtol=1e-4

f_NINO = model[:Interactions,  :f_NINO]
f_NINO_compare = benchmark.f_NINO
@test f_NINO ≈ f_NINO_compare rtol=1e-4
