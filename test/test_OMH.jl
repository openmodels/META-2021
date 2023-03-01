using Test, CSV
include("../src/basemodel.jl")
include("../src/components/OMH.jl")

benchmark = CSV.read("../data/benchmark/OMH.csv", DataFrame)

## Setup the model

model = test_model()
omh = addOMH(model, "Whiteman et al. beta 20 years")

omh[:T_AT] = benchmark.Temperature
omh[:uniforms] = 1 .- benchmark."Binomial draw"

run(model)

## Test the model

ch4 = model[:OMH,  :CH4_OMH]
ch4_compare = benchmark."CH4_OMH (MtCH4)"

@test ch4 ≈ ch4_compare rtol=1e-4
