using Test, CSV
include("../src/basemodel.jl")

benchmark = CSV.read("../data/benchmark/PCF.csv", DataFrame)

## Setup the model

model = test_model()
include("../src/components/PCF.jl")

pcfmodel = addPCFModel(model, "Fit of Hope and Schaefer (2016)")
pcfmodel[:T_AT] = benchmark."Atmospheric temperature"
# pcfmodel[:year] = benchmark."Year"
# pcfmodel[:CO2_PF] = benchmark."CO2_PF"
# pcfmodel[:CH4_PF] = benchmark."CH4_PF"

run(model)

## Test the model

CO2_PF = model[:PCFModel, :CO2_PF]
CO2_PF_compare = benchmark."CO2_PF"

@test CO2_PF ≈ CO2_PF_compare rtol=1e-4
