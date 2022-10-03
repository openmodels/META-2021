using Test, CSV, Mimi, Random, Distributions
include("../src/basemodel.jl")
include("../src/components/CH4Model.jl")
include("../src/components/RCP.jl")

benchmark = CSV.read("../data/benchmark/CH4Model.csv", DataFrame) #Generated with RCP8.5

## Set up the model

model = test_model()
RCPmodel = addRCP(model, "RCP8.5") #RCP8.5 to match benchmark file
CH4model = addCH4Model(model, "Value")

CH4model[:ch4_rcp] = RCPmodel[:ch4_rcp]
CH4model[:ch4_conc_rcp] = RCPmodel[:ch4_conc_rcp]
CH4model[:n2o_conc_rcp] = RCPmodel[:n2o_conc_rcp]

run(model)

## Test the model

concentration_CH4 = model[:CH4Model, :CH4_concentration]
concentration_CH4_compare = benchmark."Concentration"

forcing_CH4 = model[:CH4Model, :F_CH4]
forcing_CH4_compare = benchmark."Forcing (W/m2)"

@test concentration_CH4 ≈ concentration_CH4_compare rtol=1e-4
@test forcing_CH4 ≈ forcing_CH4_compare rtol=1e-4