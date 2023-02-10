using Test, CSV, DataFrames
using HypothesisTests
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")
include("../src/lib/presets.jl")

do_test = "notp"

if do_test == "full"
    model = full_model(rcp="RCP4.5", ssp="SSP2")
    update_param!(model, :CH4Model, :decay_rate, 1 / 12.4)
    prepare_montecarlo!(model)
    sim = getsim(model, "Fit of Hope and Schaefer (2016)", # PCF
                 "Cai et al. central value", # AMAZ
                 "Nordhaus central value", # GIS
                 "Distribution", # WAIS
                 "Distribution", # SAF
                 true, # ism_used
                 true, # omh_used
                 true, # amoc_used
                 false, # persit
                 false, # emuc
                 false) # prtp
elseif do_test == "some" # PCFGISISMSAF
    model = full_model(rcp="RCP4.5", ssp="SSP2"; interaction=false, omh=false, amaz=false, wais=false, amoc=false)
    update_param!(model, :CH4Model, :decay_rate, 1 / 12.4)
    prepare_montecarlo!(model)
    sim = getsim(model, "Fit of Hope and Schaefer (2016)", # PCF
                 "none", # AMAZ
                 "Nordhaus central value", # GIS
                 "none", # WAIS
                 "Distribution", # SAF
                 true, # ism_used
                 false, # omh_used
                 false, # amoc_used
                 false, # persit
                 false, # emuc
                 false) # prtp
elseif do_test == "notp"
    model = base_model(rcp="RCP4.5")
    update_param!(model, :CH4Model, :decay_rate, 1 / 12.4)
    ## prepare_montecarlo!(model)
    draws = getsim_base(500, false, false, false)
    results = runsim_base(model, draws)
end

si = run(sim, model, 500)

if do_test == "full"
    benchmark = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)
else
    benchmark = CSV.read("../data/benchmark/ExcelMETA-notp.csv", DataFrame)
end

mapping = Dict{String, Tuple{Symbol, Symbol}}("Nonlinear SAF / Distribution" => (:SAFModel, :delta),
                                              "Warming half-life / Distribution" => (:SAFModel, :FRT),
                                              "xi_1" => (:TemperatureModel, :xi_1),
                                              "F_2XCO2" => (:TemperatureModel, :F_2xCO2),
                                              "T_2xCO2" => (:TemperatureModel, :fair_ECS),
                                              "gamma" => (:TemperatureModel, :fair_gamma),
                                              "C0" => (:TemperatureModel, :fair_C_0),
                                              "Dataset 1" => (:CH4Model, :ch4_alpha))

tvals = Dict{String, Float64}()
missings = String[]
errors = String[]
for col in names(benchmark)[2:end]
    println(col)
    info = split(col, " / ")
    try
        if col ∈ ["Global welfare", "World consumption per capita in 2020"]
            continue
        elseif col ∈ ["a_0", "a_1", "a_3", "rho_1", "rho_2", "rho_3"]
            df = getdataframe(si, :CO2Model, Symbol(replace(col, "_" => "")))
        elseif col ∈ ["r_{0) / Distribution", "r_{C} / Distribution", "r_{T} / Distribution"]
            df = getdataframe(si, :PostTemperature, Symbol(replace(info[1], "{" => "", "}" => "", ")" => "")))
        elseif col ∈ keys(mapping)
            df = getdataframe(si, mapping[col][1], mapping[col][2])
        elseif info[2] == "Kessler probabilistic" || info[1] ∈ ["Delta_AMAZ", "avoldot0"]
            continue # not included in test
        elseif info[2] == "Atmospheric temperature"
            df = getdataframe(si, :TemperatureModel, :T_AT)
        elseif info[2] == "SLR total (m)"
            df = getdataframe(si, :SLRModel, :SLR)
        elseif info[2] == "World aggregate consumption per capita"
            continue
            # df = getdataframe(si, :Utility, :equiv_conspc)
        elseif info[2] ∈ ["beta1dist", "beta2dist", "errdist"]
            continue # not saved
        elseif info[1] == "waisrate"
            df = getdataframe(si, :WAISmodel, :waisrate)
            df = DataFrame(:waisrate => df.waisrate * 1000, :trialnum => df.trialnum)
        else
            push!(missings, col)
            continue
        end
        if length(info) == 2 && info[2] != "Distribution"
            vals = convert(Vector{Float64}, df[df.time .== parse(Int64, info[1]), 2])
        else
            vals = convert(Vector{Float64}, df[!, 1])
        end
        comps = benchmark[!, col]

        result = UnequalVarianceTTest(vals[.!isnan.(vals)], comps)
        tvals[col] = result.t
    catch
        push!(errors, col)
    end
    # @test result.t < 2.56 # 99%
end
