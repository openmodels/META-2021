using Test, CSV, Random, Distributions
include("../src/basemodel.jl")

benchmark_df = CSV.read("../data/benchmark/ISM.csv", DataFrame)

## Setup the model

model = test_model()
include("../src/components/ISM.jl")

ismmodel = addISMModel(model, "Value")
ismmodel[:st_ppm] = benchmark_df."S_t (ppm)"
ismmodel[:SO_2] = benchmark_df."Emissions|Sulfur"
ismmodel[:T_AT] = benchmark_df."Atmospheric temperature"

benchmark = Matrix(benchmark_df)
ismmodel[:uniforms] = benchmark[1:end, 6:end]

run(model)

## Test the model

D_ISM = model[:ISMModel, :D_ISM]
D_ISM_compare = benchmark_df."Damages"

@test D_ISM ≈ D_ISM_compare rtol=1e-4
