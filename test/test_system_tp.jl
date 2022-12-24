using Test, CSV, DataFrames
include("../src/MimiMETA.jl")
include("../src/lib/presets.jl")

benchmark = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)

## Setup the model

alltmaxdiff = []
allsmaxdiff = []
allgweldiff = []
for rr in 1:nrow(benchmark)
    println(rr)
    global model = full_model() # XXX: Consumption currently refers to this...
    preset_fill_tp(model, benchmark, rr)

    run(model)

    ## Test the model

    T_AT = model[:TemperatureModel, :T_AT][11:10:191]
    T_AT_compare = collect(benchmark[rr, 2:20])

    @test maximum(abs.(T_AT .- T_AT_compare)) < .1

    SLR = model[:SLRModel, :SLR][11:10:191]
    SLR_compare = collect(benchmark[rr, 21:39])

    @test SLR ≈ SLR_compare atol=1e-1

    globalwelfare = sum(model[:Utility, :world_disc_utility][11:191])
    globalwelfare_compare = benchmark."Global welfare"[rr]

    @test globalwelfare ≈ globalwelfare_compare rtol=1e-2

    push!(alltmaxdiff, (T_AT .- T_AT_compare)[findmax(abs.(T_AT .- T_AT_compare))[2]])
    push!(allsmaxdiff, (SLR .- SLR_compare)[findmax(abs.(SLR .- SLR_compare))[2]])
    push!(allgweldiff, globalwelfare - globalwelfare_compare)
end

df = DataFrame(:tmaxdiff => alltmaxdiff, :smaxdiff => allsmaxdiff, :gweldiff => allgweldiff)
CSV.write("errors-tp.csv", df)
