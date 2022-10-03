using Test, CSV
include("../src/basemodel.jl")

benchmark = CSV.read("../data/benchmark/PatternScaling.csv", DataFrame)

## Setup the model

include("../src/components/PatternScaling.jl")

model = test_model()
pattscale = addPatternScaling(model)

pattscale[:T_AT] = benchmark."Atmospheric temperature" # Use the final (post-SAF) temperatures

run(model)

## Test the model

scale = model[:PatternScaling,  :T_country]
scale_compare = Matrix(benchmark[!, 3:end-1])

@test scale ≈ scale_compare
