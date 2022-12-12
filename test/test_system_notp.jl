using Test, CSV, DataFrames
include("../src/MimiMETA.jl")
include("../src/lib/presets.jl")

benchmark = CSV.read("../data/benchmark/ExcelMETA-notp.csv", DataFrame)

## Setup the model
model = base_model(; rcp="RCP4.5")

for rr in 1:nrow(benchmark)
    println(rr)
    preset_fill_notp(model, benchmark, rr)

    run(model)

    ## Test the model

    T_AT = model[:TemperatureModel, :T_AT][11:10:191]
    T_AT_compare = collect(benchmark[rr, 2:20])

    @test T_AT ≈ T_AT_compare atol=1e-4

    SLR = model[:SLRModel, :SLR][11:10:191]
    SLR_compare = collect(benchmark[rr, 21:39])

    @test SLR ≈ SLR_compare atol=1e-4

    globalwelfare = sum(model[:Utility, :world_disc_utility][11:191])
    globalwelfare_compare = benchmark."Global welfare"[rr]

    @test globalwelfare ≈ globalwelfare_compare rtol=1e-4
end
