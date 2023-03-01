using Test, CSV
include("../src/basemodel.jl")

benchmark = CSV.read("../data/benchmark/SAF.csv", DataFrame)

## Setup the model

model = test_model()
include("../src/components/SAF.jl")

safmodel = addSAFModel(model, "Distribution mean")
safmodel[:F] = benchmark."Forcing"

run(model)

## Test the model

T_AT_adjustment = model[:SAFModel, :T_AT_adjustment]
T_AT_adjustment_compare = benchmark."Temp Correction"

@test T_AT_adjustment ≈ T_AT_adjustment_compare rtol=1e-4
