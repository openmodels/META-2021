using Test, CSV
include("../src/basemodel.jl")
include("../src/components/NonMarketDamages.jl")

benchmark = CSV.read("../data/benchmark/NonMarket-loss.csv", DataFrame)
benchmark_temp = CSV.read("../data/benchmark/TemperatureModel-May2022.csv", DataFrame)
benchmark_cons = CSV.read("../data/benchmark/Consumption.csv", DataFrame)

## Setup the model

model = test_model()
nonmarket = addNonMarketDamages(model)

nonmarket[:T_AT] = benchmark_temp."Atmospheric temperature"
nonmarket[:conspc] = Matrix(benchmark_cons[2:end, 2:end-5])

run(model)

## Test the model

lossfactor = convert(Matrix{Float64}, model[:NonMarketDamages,  :lossfactor])
lossfactor_compare = Matrix(benchmark[!, 2:end])

@test lossfactor ≈ lossfactor_compare rtol=1e-4
