using Test, CSV
include("../src/basemodel.jl")

benchmark_cons = CSV.read("../data/benchmark/Consumption.csv", DataFrame)
benchmark_lossfactor = CSV.read("../data/benchmark/NonMarket-loss.csv", DataFrame)

benchmark = CSV.read("../data/benchmark/Utility.csv", DataFrame)

## Setup the model

model = test_model()
include("../src/components/Utility.jl")

utility = addUtility(model, "SSP2")

utility[:conspc] = Matrix(benchmark_cons[2:192, 2:195])
utility[:lossfactor] = Matrix(benchmark_lossfactor[!, 2:end])

run(model)

## Test the model

world_disc_utility = model[:Utility, :world_disc_utility]
world_disc_utility_compare = benchmark."World discounted utility"

@test world_disc_utility ≈ world_disc_utility_compare rtol=1e-4
