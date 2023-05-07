using Test, CSV, Mimi, Random, Distributions

include("../src/basemodel.jl")
include("../src/components/AIS.jl")

## benchmark = CSV.read("../data/benchmark/AIS.csv", DataFrame)

## Set up the model

model = test_model(startyear=1750)
AISmodel = addAISmodel(model)

#Add temperature time series from benchmark file

temperature_data = CSV.read("../data/benchmark/TemperatureModel-May2022.csv", DataFrame)
alltat = [collect(range(0, 0.854; length=2010 - 1750 + 1))[1:end-1];
          temperature_data."Atmospheric temperature_1"]

AISmodel[:T_AT] = alltat
AISmodel[:T_AT_tminus100] = [zeros(100); alltat[1:end-100]]

run(model)

## Test the model

SLR_AIS = model[:AISmodel, :SLR_AIS]
SLR_AIS_compare = benchmark."Cumulative Dynamic and SMB SLR"

@test SLR_AIS ≈ SLR_AIS_compare rtol=1e-4
