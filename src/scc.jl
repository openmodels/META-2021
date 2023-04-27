using Mimi
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")

function calculate_scc(model::Model, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    mm = calculate_scc_setup(model, pulse_year, pulse_size)
    run(mm)
    calculate_scc_marginal(mm, pulse_year, emuc)
end

function calculate_scc_setup(model::Model, pulse_year::Int64, pulse_size::Float64)
    mm = create_marginal_model(model, pulse_size)

    pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)
    run(mm)

    mm.modified[:CO2Converter, :co2_extra][pulse_index] = pulse_size

    mm
end

function calculate_scc_marginal(mm::Union{MarginalModel, MarginalInstance}, pulse_year::Int64, emuc::Float64)
    pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)

    globalwelfare_marginal = sum(mm[:Utility, :world_disc_utility][pulse_index:dim_count(model, :time)])

    global_conspc = sum(mm.base[:Consumption, :conspc][pulse_index, :] .* mm.base[:Utility, :pop][pulse_index, :]) / mm.base[:Utility, :world_population][pulse_index]
    -(globalwelfare_marginal / (global_conspc^-emuc)) / 1e9
end

function calculate_scc_base_mc(model::Model, trials::Int64, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    mm = calculate_scc_setup(model, pulse_year, pulse_size)

    function setsim_base_scc(inst::Union{ModelInstance, MarginalInstance}, draws::DataFrame, ii::Int64, pulse_year::Int64, pulse_size::Float64)
        setsim_base(inst, draws, ii)
        pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)
        inst.modified[:CO2Converter, :co2_extra][pulse_index] = pulse_size
    end

    sim_base(mm, trials, persist_dist, emuc_dist, prtp_dist; save_rvs=false,
             setsim=(inst, draws, ii) -> setsim_base_scc(inst, draws, ii, pulse_year, pulse_size),
             getsim=(inst, draws; save_rvs) -> calculate_scc_marginal(inst, pulse_year, emuc))
end

# model = base_model(; rcp="CP-Base", tdamage="pointestimate", slrdamage="mode")
# calculate_scc(model, 2020, 10., 1.5)
# sccs = calculate_scc_base_mc(model, 100, false, false, false, 2020, 10., 1.5)
# [mean(sccs[:other]), std(sccs[:other]), median(sccs[:other])]

function calculate_scc_full_mc(model::Model, trials::Int64, pcf_calib::String, amazon_calib::String, gis_calib::String, wais_calib::String, saf_calib::String, ais_dist::Bool, ism_used::Bool, omh_used::Bool, amoc_used::Bool, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    mm = calculate_scc_setup(model, pulse_year, pulse_size)

    function setsim_full_scc(inst::Union{ModelInstance, MarginalInstance}, draws::DataFrame, ii::Int64, ism_used::Bool, omh_used::Bool, amoc_used::Bool, amazon_calib::String, wais_calib::String, ais_dist::Bool, pulse_year::Int64, pulse_size::Float64)
        setsim_full(inst, draws, ii, ism_used, omh_used, amoc_used, amazon_calib, wais_calib, ais_dist)
        pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)
        inst.modified[:CO2Converter, :co2_extra][pulse_index] = pulse_size
    end

    sim_full(mm, trials, pcf_calib, amazon_calib, gis_calib, wais_calib, saf_calib,
             ais_dist, ism_used, omh_used, amoc_used, persist_dist, emuc_dist, prtp_dist; save_rvs=false,
             setsim=(inst, draws, ii, ism_used, omh_used, amoc_used, amazon_calib, wais_calib, ais_dist) -> setsim_full_scc(inst, draws, ii, ism_used, omh_used, amoc_used, amazon_calib, wais_calib, ais_dist, pulse_year, pulse_size),
             getsim=(inst, draws; save_rvs) -> calculate_scc_marginal(inst, pulse_year, emuc))
end

# model = full_model(rcp="RCP4.5", ssp="SSP2")
# calculate_scc(model, 2020, 1.0, 1.5)
# sccs = calculate_scc_full_mc(model, 100,
#                              "Fit of Hope and Schaefer (2016)", # PCF
#                              "Cai et al. central value", # AMAZ
#                              "Nordhaus central value", # GIS
#                              "Distribution", # WAIS
#                              "Distribution", # SAF
#                              false, # ais_used
#                              true, # ism_used
#                              true, # omh_used
#                              true, # amoc_used
#                              false, # persit
#                              false, # emuc
#                              false, # prtp
#                              2020, 10., 1.5)
# [mean(sccs[:other]), std(sccs[:other]), median(sccs[:other])]
