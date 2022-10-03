using Test, CSV
include("../src/basemodel.jl")

benchmark = CSV.read("../data/benchmark/TemperatureModel.csv", DataFrame)

## Setup the model

model = test_model()
include("../src/components/TemperatureModel.jl")
include("../src/components/Forcing.jl")

forcing = addForcing(model, "Best fit multi-model mean")
forcing[:st_ppm] = benchmark."S_t (ppm)"
forcing[:F_CH4] = benchmark."Forcing|CH4"
forcing[:F_EX] = benchmark."Forcing|other agents"

temperaturemodel = addTemperatureModel(model, "Best fit multi-model mean")

temperaturemodel[:F] = forcing[:F]

run(model)

## Test the model

T_AT = model[:TemperatureModel, :T_AT]
T_AT_compare = benchmark."Atmospheric temperature"

@test T_AT ≈ T_AT_compare rtol=1e-4
