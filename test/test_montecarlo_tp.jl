using Test, CSV, DataFrames
using HypothesisTests
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")
include("../src/lib/presets.jl")

mapping = Dict{String, Tuple{Symbol, Symbol}}("Nonlinear SAF / Distribution" => (:SAFModel, :saf_delta),
                                              "Warming half-life / Distribution" => (:SAFModel, :FRT),
                                              "xi_1" => (:TemperatureModel, :xi_1),
                                              "F_2XCO2" => (:TemperatureModel, :F_2xCO2),
                                              "T_2xCO2" => (:TemperatureModel, :fair_ECS),
                                              "gamma" => (:TemperatureModel, :fair_gamma),
                                              "C0" => (:TemperatureModel, :fair_C_0),
                                              "Dataset 1" => (:CH4Model, :ch4_alpha))

for do_test in ["notp", "full", "some"]
    if do_test == "full"
        ## Run a test with all TPs
        global model = full_model(rcp="RCP4.5", ssp="SSP2")
        draws = getsim(500, "Fit of Hope and Schaefer (2016)", # PCF
                       "Cai et al. central value", # AMAZ
                       "Nordhaus central value", # GIS
                       "Distribution", # WAIS
                       "Distribution", # SAF
                       false, # persit
                       false, # emuc
                       false) # prtp
        results = runsim(model, draws, true, # ism_used
                         true, # omh_used
                         true, # amoc_used
                         "Cai et al. central value", # AMAZ
                         "Distribution") # WAIS
    elseif do_test == "some" # PCFGISISMSAF
        ## Run a test with some TPs
        global model = full_model(rcp="RCP4.5", ssp="SSP2"; interaction=false, omh=false, amaz=false, wais=false, amoc=false)
        draws = getsim(500, "Fit of Hope and Schaefer (2016)", # PCF
                       "none", # AMAZ
                       "Nordhaus central value", # GIS
                       "none", # WAIS
                       "Distribution", # SAF
                       false, # persit
                       false, # emuc
                       false) # prtp
        results = runsim(model, draws, true, # ism_used
                         false, # omh_used
                         false, # amoc_used
                         "none", # AMAZ
                         "none") # WAIS
    elseif do_test == "notp"
        ## Run a test with no TPs
        global model = base_model(rcp="RCP4.5")
        draws = getsim_base(500, false, false, false)
        results = runsim_base(model, draws)
    end

    run(model) # run once for simdataframe

    for compare_same in [true, false]
        if (compare_same && do_test == "full") || (!compare_same && do_test == "some")
            benchmark = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)
        elseif (compare_same && do_test == "some") || (!compare_same && do_test == "notp")
            benchmark = CSV.read("../data/benchmark/ExcelMETA-PCFGISISMSAF.csv", DataFrame)
        else
            benchmark = CSV.read("../data/benchmark/ExcelMETA-notp.csv", DataFrame)
        end

        tvals = Dict{String, Float64}()
        missings = String[]
        errors = String[]
        for col in names(benchmark)[2:end]
            info = split(col, " / ")
            try
                if col ∈ ["Global welfare", "World consumption per capita in 2020"]
                    continue
                elseif col ∈ ["a_0", "a_1", "a_3", "rho_1", "rho_2", "rho_3"]
                    df = simdataframe(model, results, :CO2Model, Symbol(replace(col, "_" => "")))
                elseif col ∈ ["r_{0) / Distribution", "r_{C} / Distribution", "r_{T} / Distribution"]
                    df = simdataframe(model, results, :PostTemperature, Symbol(replace(info[1], "{" => "", "}" => "", ")" => "")))
                elseif col ∈ keys(mapping)
                    df = simdataframe(model, results, mapping[col][1], mapping[col][2])
                elseif info[2] == "Kessler probabilistic" || info[1] ∈ ["Delta_AMAZ", "avoldot0"]
                    continue # not included in test
                elseif info[2] == "Atmospheric temperature"
                    df = simdataframe(model, results, :TemperatureModel, :T_AT)
                elseif info[2] == "SLR total (m)"
                    df = simdataframe(model, results, :SLRModel, :SLR)
                elseif info[2] == "World aggregate consumption per capita"
                    continue
                    # df = simdataframe(model, results, :Utility, :equiv_conspc)
                elseif info[2] ∈ ["beta1dist", "beta2dist", "errdist"]
                    continue # not saved
                elseif info[1] == "waisrate"
                    df = simdataframe(model, results, :WAISmodel, :waisrate)
                    df = DataFrame(:trialnum => df.trialnum, :waisrate => df.waisrate * 1000)
                else
                    push!(missings, col)
                    continue
                end
                if length(info) == 2 && info[2] != "Distribution"
                    vals = convert(Vector{Float64}, df[df.time .== parse(Int64, info[1]), 2])
                else
                    vals = convert(Vector{Float64}, df[!, 2])
                end
                comps = benchmark[!, col]

                result = UnequalVarianceTTest(vals[.!isnan.(vals)], comps)
                tvals[col] = result.t
                if compare_same
                    @test abs(result.t) < 2.56 # 99%
                end
            catch
                push!(errors, col)
            end
        end

        if !compare_same
            @test maximum(abs.(values(tvals))) > 2.56
        end
    end
end
