using Test, CSV
include("../src/basemodel.jl")

benchmark = CSV.read("../data/benchmark/CO2Model.csv", DataFrame)

## Setup the model

model = test_model()
include("../src/components/CO2Model.jl")

co2model = addCO2Model(model, "Expectation")
co2model[:co2_rcp] = [benchmark."Emissions (GtCO2)"[2:end]; 0] # last value unused
co2model[:alpha] = benchmark."alpha"

run(model)

## Test the model

st = model[:CO2Model,  :st_ppm]
st_compare = benchmark."S_t (ppm)"

@test st ≈ st_compare rtol=1e-4
