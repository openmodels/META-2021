using Test, CSV
include("../src/basemodel.jl")
include("../src/components/AMOC.jl")

benchmark = CSV.read("../data/benchmark/Hosing.csv", DataFrame)
benchmark_scale = CSV.read("../data/benchmark/PatternScaling_AMOConly.csv", DataFrame)

## Setup the model

model = test_model()
amoc = addAMOC(model, "Hadley")

amoc[:T_AT] = benchmark.Temperature[2:end]
amoc[:f_AMOC] = repeat([1], dim_count(model, :time))
amoc[:uniforms] = 1 .- benchmark."Binomial draw"[2:end]

amoc[:scale_country] = Matrix(benchmark_scale[!, 2:end-1])

run(model)

## Test the model

deltaT = model[:AMOC, :deltaT_country_AMOC]
deltaT_compare = Matrix(benchmark[2:end, 6:end])

@test deltaT[:, 1:34] ≈ deltaT_compare[:, 1:34]

benchmark_temps = CSV.read("../data/benchmark/AMOC_Temperature.csv", DataFrame)

temps = model[:AMOC,  :T_country_AMOC]
temps_compare = Matrix(benchmark_temps[!, 2:end])

@test temps[:, 1:34] ≈ temps_compare[:, 1:34]
