using Test, CSV
include("../src/basemodel.jl")
include("../src/components/Consumption.jl")

benchmark_temp = CSV.read("../data/benchmark/AMOC_Temperature-May2022.csv", DataFrame)
benchmark_slr = CSV.read("../data/benchmark/SLRModel-May2022.csv", DataFrame)

## Setup the model
include("../src/components/Consumption.jl")

model = test_model()
cons = addConsumption(model, "none", "none", "SSP2")

cons[:T_country] = Matrix(benchmark_temp[!, 2:end])
cons[:SLR] = benchmark_slr[!, 4]

run(model)

## Test the model without impacts

gdppc = model[:Consumption, :gdppc]
gdppc_compare = Matrix(CSV.read("../data/benchmark/National gdp per capita ppp.csv", DataFrame, limit=290)[2:end, 2:end])

@test gdppc[1:3, 1:3] ≈ gdppc_compare[1:3, 1:3]
@test gdppc[1:100, 1:3] ≈ gdppc_compare[1:100, 1:3] rtol=1e-4

benchmark = CSV.read("../data/benchmark/Consumption-noimp.csv", DataFrame)

conspc = model[:Consumption, :conspc]
conspc_compare = Matrix(benchmark[2:end, 2:end])

@test conspc[1:3, 1:3] ≈ conspc_compare[1:3, 1:3]
@test conspc[1:100, 1:194] ≈ conspc_compare[1:100, 1:194] rtol=1e-4

## Test the model with impacts

model = test_model()
cons = addConsumption(model, "pointestimate", "mode", "SSP2")

cons[:T_country] = Matrix(benchmark_temp[!, 2:end])
cons[:SLR] = benchmark_slr[!, 4]

run(model)

benchmark = CSV.read("../data/benchmark/Consumption.csv", DataFrame)

conspc = model[:Consumption, :conspc]
conspc_compare = Matrix(benchmark[2:end, 2:end])

@test conspc[1:3, 1:3] ≈ conspc_compare[1:3, 1:3]
@test conspc[1:100, 1:194] ≈ conspc_compare[1:100, 1:194] rtol=1e-4
