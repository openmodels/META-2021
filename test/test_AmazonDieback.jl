using Test, CSV
include("../src/basemodel.jl")
include("../src/components/AmazonDieback.jl")

benchmark = CSV.read("../data/benchmark/AmazonDieback.csv", DataFrame)

## Setup the model

model = test_model()
amaz = addAmazonDieback(model, "Cai et al. central value")

amaz[:T_AT] = benchmark.Temperature
amaz[:uniforms] = 1 .- benchmark."Binomial draw"
amaz[:probmult] = ones(nrow(benchmark))

run(model)

## Test the model

co2 = model[:AmazonDieback,  :CO2_AMAZ]
co2_compare = benchmark."CO2_AMAZ (GtCO2)"

@test co2 ≈ co2_compare rtol=1e-4
