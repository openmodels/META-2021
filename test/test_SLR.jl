using Test, CSV, Mimi, Random, Distributions
include("../src/basemodel.jl")
include("../src/components/SLR.jl")
include("../src/components/RCP.jl")

benchmark = CSV.read("../data/benchmark/SLRModel-May2022.csv", DataFrame) 

## Set up the model

model = test_model()

SLRmodel = addSLR(model)

## Add temperature time series from benchmark file
SLRmodel[:T_AT] = benchmark."Temperature" 

run(model)

## Test the model

SLR_therm = model[:SLRModel, :SLR_therm]
SLR_therm_compare = benchmark."SLR from thermal expansion and melt from glaciers and small ice caps (m)"

@test SLR_therm ≈ SLR_therm_compare rtol=1e-4