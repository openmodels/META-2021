using Test, CSV, DataFrames
include("../src/MimiMETA.jl")
include("../src/lib/presets.jl")

benchmark = CSV.read("../data/benchmark/ExcelMETA-PCFGISISMSAF.csv", DataFrame)

## Setup the model

alltmaxdiff = []
allsmaxdiff = []
allgweldiff = []
for rr in 1:nrow(benchmark)
    println(rr)
    global model = full_model(; interaction=false, omh=false, amaz=false, wais=false, amoc=false) # XXX: Consumption currently refers to this...

    preset_fill_notp(model, benchmark, rr)

    # Not using distribution
    # myupdate_param!(model, :PCFModel, :propCH4, benchmark."propCH4 / Kessler probabilistic"[rr])
    # myupdate_param!(model, :PCFModel, :beta_PF, benchmark."beta_PF / Kessler probabilistic"[rr])
    # myupdate_param!(model, :PCFModel, :C_PF, benchmark."C_PF (GtC) / Kessler probabilistic"[rr])
    # myupdate_param!(model, :PCFModel, :propPassive, benchmark."propPassive / Kessler probabilistic"[rr])
    # myupdate_param!(model, :PCFModel, :tau, benchmark."tau (years) / Kessler probabilistic"[rr])

    # Not using distribution
    # myupdate_param!(model, :GISModel, :avoldot0, benchmark."avoldot0 / Distribution"[rr])

    myupdate_param!(model, :SAFModel, :ECS, benchmark.T_2xCO2[rr])
    if benchmark."Nonlinear SAF / Distribution"[rr] != "Error"
        myupdate_param!(model, :SAFModel, :saf_delta, benchmark."Nonlinear SAF / Distribution"[rr])
    end
    if benchmark."Warming half-life / Distribution"[rr] != "Error"
        myupdate_param!(model, :SAFModel, :FRT, benchmark."Warming half-life / Distribution"[rr])
    end

    raindrawstart = findfirst(names(benchmark) .== "2010 / Day")
    bindrawstarts = findall(x -> occursin("2010 / Binomial draw", x), names(benchmark))
    raindrawend = bindrawstarts[4] - 1
    draws = reshape(collect(benchmark[rr, raindrawstart:raindrawend]), (dim_count(model, :monsoonsteps), dim_count(model, :time)))'
    myupdate_param!(model, :ISMModel, :uniforms, draws)

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
CSV.write("errors-PCFGISISMSAF.csv", df)
