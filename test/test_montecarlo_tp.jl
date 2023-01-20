using Test, CSV, DataFrames
using HypothesisTests
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")
include("../src/lib/presets.jl")

model = full_model(rcp="RCP4.5", ssp="SSP2")
update_param!(model, :CH4Model, :decay_rate, 1 / 12.4)
prepare_montecarlo!(model)
sim = getsim("Fit of Hope and Schaefer (2016)", "Cai et al. central value", "Nordhaus central value", "Distribution", "Distribution", false, false, false)

si = run(sim, model, 500)

benchmark = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)

tvals = Float64[]
for col in names(benchmark)[2:end]
    println(col)
    info = split(col, " / ")
    if length(info) == 2
        if info[2] == "Atmospheric temperature"
            df = getdataframe(si, :TemperatureModel, :T_AT)
        elseif info[2] == "SLR total (m)"
            df = getdataframe(si, :SLRModel, :SLR)
        elseif info[2] == "World aggregate consumption per capita"
            continue
            # df = getdataframe(si, :Utility, :equiv_conspc)
        else
            break
        end
        vals = convert(Vector{Float64}, df[df.time .== parse(Int64, info[1]), 2])
        comps = benchmark[!, col]

        result = UnequalVarianceTTest(vals, comps)
        push!(tvals, result.t)
        # @test result.t < 2.56 # 99%
    else
        if col ∈ ["Global welfare", "World consumption per capita in 2020"]
            continue
        else
            break
        end
    end
end
