using Mimi
using Distributions
using CSV
using DataFrames
using LinearAlgebra
using ProgressMeter
using MimiFAIRv2
include("../src/lib/presets.jl")
include("../src/lib/MimiFAIR_monte_carlo.jl")

import Mimi.ModelInstance

function make_lognormal(riskmu, risksd)
    mu = log(riskmu^2 / sqrt(risksd^2 + riskmu^2))
    sd = sqrt(log(1 + (risksd /  riskmu)^2))
    LogNormal(mu, sd)
end

# Master function for base model (uses helpers below)
function sim_base(model::Model, trials::Int64, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool; save_rvs::Bool=true)
    draws = presim_base(trials, persist_dist, emuc_dist, prtp_dist)

    sim = create_fair_monte_carlo(model, trials; end_year=2200,
                                  data_dir=joinpath(dirname(pathof(MimiFAIRv2)), "..", "data",
                                                    "large_constrained_parameter_files"),
                                  delete_downloaded_data=false,
                                  other_mc_set=(inst, ii) -> setsim_base(inst, draws, ii),
                                  other_mc_get=(inst) -> getsim_base(inst, draws, save_rvs=save_rvs))
    sim()
end


function presim_base(trials::Int64, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool)
    draws = DataFrame(mc=1:trials)

    # Utility

    if persist_dist
        draws.Consumption_damagepersist = rand(Uniform(0, 1), trials)
    end
    if emuc_dist
        draws.Consumption_EMUC = rand(TriangularDist(0.5, 2, 1.5), trials)
    end
    if prtp_dist
        draws.Consumption_PRTP = rand(TriangularDist(0.001, 0.02, 0.01), trials)
    end

    draws
end

function setsim_base(inst::ModelInstance, draws::DataFrame, ii::Int64)
    for jj in 2:ncol(draws)
        update_param!(inst, Symbol(names(draws)[jj]), draws[ii, jj])
    end

    # Damages

    update_param!(inst, :Consumption_seeds, rand(DiscreteUniform(1, typemax(Int64)), dim_count(model, :country)))
    update_param!(inst, :Consumption_slruniforms, rand(Uniform(0, 1), dim_count(model, :country)))
end

function getsim_base(inst::ModelInstance, draws::DataFrame; save_rvs::Bool=true)
    mcres = Dict{Symbol, Any}()
    mcres[:temperature_T] = copy(inst[:temperature, :T])
    mcres[:SLRModel_SLR] = copy(inst[:SLRModel, :SLR])
    mcres[:Utility_world_disc_utility] = inst[:Utility, :world_disc_utility]

    if save_rvs
        for jj in 2:ncol(draws)
            mcres[Symbol(names(draws)[jj])] = draws[ii, jj]
        end
    end

    mcres
end

function presim(trials::Int64, pcf_calib::String, amazon_calib::String, gis_calib::String, wais_calib::String, saf_calib::String, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool)
    draws = presim_base(trials, persist_dist, emuc_dist, prtp_dist)

    # Interactions

    dists = allinteractrates((symbol, ratemu, ratese) -> Normal(ratemu, ratese))

    draws.Interactions_gis2amoc = rand(dists["gis"]["amoc"], trials)
    draws.Interactions_wais2amoc = rand(dists["wais"]["amoc"], trials)
    draws.Interactions_amaz2amoc = rand(dists["amaz"]["amoc"], trials)
    draws.Interactions_nino2amoc = rand(dists["nino"]["amoc"], trials)
    draws.Interactions_amoc2gis = rand(dists["amoc"]["gis"], trials)
    draws.Interactions_wais2gis = rand(dists["wais"]["gis"], trials)
    draws.Interactions_amaz2gis = rand(dists["amaz"]["gis"], trials)
    draws.Interactions_nino2gis = rand(dists["nino"]["gis"], trials)
    draws.Interactions_amoc2wais = rand(dists["amoc"]["wais"], trials)
    draws.Interactions_gis2wais = rand(dists["gis"]["wais"], trials)
    draws.Interactions_amaz2wais = rand(dists["amaz"]["wais"], trials)
    draws.Interactions_nino2wais = rand(dists["nino"]["wais"], trials)
    draws.Interactions_amoc2amaz = rand(dists["amoc"]["amaz"], trials)
    draws.Interactions_gis2amaz = rand(dists["gis"]["amaz"], trials)
    draws.Interactions_wais2amaz = rand(dists["wais"]["amaz"], trials)
    draws.Interactions_nino2amaz = rand(dists["nino"]["amaz"], trials)
    draws.Interactions_amoc2nino = rand(dists["amoc"]["nino"], trials)
    draws.Interactions_gis2nino = rand(dists["gis"]["nino"], trials)
    draws.Interactions_wais2nino = rand(dists["wais"]["nino"], trials)
    draws.Interactions_amaz2nino = rand(dists["amaz"]["nino"], trials)

    # AIS

    draws.AIS_ω = rand(Normal(-0.05, 0.004), trials)
    draws.AIS_λ = rand(Uniform(7, 16), trials)

    # Permafrost

    if pcf_calib == "Kessler probabilistic"
        draws.PCFModel_propCH4 = rand(Normal(0.023, 0.006), trials)
        draws.PCFModel_beta_PF = rand(Normal(0.172, 0.0261), trials)
        draws.PCFModel_C_PF = rand(Normal(1035, 76.53), trials)
        draws.PCFModel_propPassive = rand(Normal(0.4, 0.055), trials)
        draws.PCFModel_tau = rand(Normal(70, 30), trials)
    end

    # Amazon

    if amazon_calib == "Distribution"
        draws.AmazonDieback_Delta_AMAZ = rand(TriangularDist(10, 250, 50), trials)
    end

    # GIS

    if gis_calib == "Distribution"
        draws.GISModel_avoldot = rand(Normal(-0.0000106, 0.0000244/0.5/100), trials) # Only works with a meltmult of 1
    end

    # WAIS

    if wais_calib == "Distribution"
        draws.WAISmodel_waisrate = rand(make_lognormal(3.3 / 1000, 1.65 / 1000), trials)
    end

    # SAF

    if saf_calib == "Distribution"
        draws.SAFModel_saf_delta = rand(TriangularDist(-1, 1, 0), trials)
        draws.SAFModel_FRT = rand(TriangularDist(10, 55, 20), trials)
    end

    draws
end

function runsim(model::Model, draws::DataFrame, ism_used::Bool, omh_used::Bool, amoc_used::Bool, amazon_calib::String, wais_calib::String; save_rvs::Bool=true)
    if wais_calib == "Distribution"
        set_param!(model, :WAISmodel, :waisrate, :WAISmodel_waisrate, 0.0033) # set up global connection
    end

    # AIS
    aisgcms = CSV.read("../data/Basal_melt_models.csv", DataFrame)
    aisresponse_EAIS = CSV.read("../data/Response functions - EAIS.csv", DataFrame)
    aisresponse_Ross = CSV.read("../data/Response functions - Ross.csv", DataFrame)
    aisresponse_Amundsen = CSV.read("../data/Response functions - Amundsen.csv", DataFrame)
    aisresponse_Weddell = CSV.read("../data/Response functions - Weddell.csv", DataFrame)
    aisresponse_Peninsula = CSV.read("../data/Response functions - Peninsula.csv", DataFrame)

    inst = Mimi.build(model)

    progress = Progress(nrow(draws), 1)

    results = Dict{Symbol, Any}[]
    for ii in 1:nrow(draws)
        for jj in 2:ncol(draws)
            if has_parameter(model.md, Symbol(names(draws)[jj]))
                update_param!(inst, Symbol(names(draws)[jj]), draws[ii, jj])
            end
        end

        # AIS

        gcmchoice = rand(DiscreteUniform(1, 19), 1)
        update_param!(inst, :AIS_β_EAIS, aisgcms.EAIS_beta[gcmchoice])
        update_param!(inst, :AIS_δ_EAIS, aisgcms.EAIS_delta[gcmchoice])
        update_param!(inst, :AIS_β_EAIS, aisgcms.Ross_beta[gcmchoice])
        update_param!(inst, :AIS_δ_EAIS, aisgcms.Ross_delta[gcmchoice])
        update_param!(inst, :AIS_β_EAIS, aisgcms.Amundsen_beta[gcmchoice])
        update_param!(inst, :AIS_δ_EAIS, aisgcms.Amundsen_delta[gcmchoice])
        update_param!(inst, :AIS_β_EAIS, aisgcms.Weddell_beta[gcmchoice])
        update_param!(inst, :AIS_δ_EAIS, aisgcms.Weddell_delta[gcmchoice])
        update_param!(inst, :AIS_β_EAIS, aisgcms.Peninsula_beta[gcmchoice])
        update_param!(inst, :AIS_δ_EAIS, aisgcms.Peninsula_delta[gcmchoice])
        icechoice = rand(DiscreteUniform(1, 17), 1)
        update_param!(inst, :AIS_R_functions_EAIS, aisresponse_EAIS[!, icechoice])
        update_param!(inst, :AIS_R_functions_Ross, aisresponse_Ross[!, icechoice])
        update_param!(inst, :AIS_R_functions_Amundsen, aisresponse_Amundsen[!, icechoice])
        update_param!(inst, :AIS_R_functions_Weddell, aisresponse_Weddell[!, icechoice])
        update_param!(inst, :AIS_R_functions_Peninsula, aisresponse_Peninsula[!, icechoice])

        # SAF

        update_param!(inst, :SAFModel_ECS, draws.TemperatureModel_fair_ECS[ii]) # SAFModel.ECS

        # Damages

        update_param!(inst, :Consumption_seeds, rand(DiscreteUniform(1, typemax(Int64)), dim_count(model, :country)))
        update_param!(inst, :Consumption_slruniforms, rand(Uniform(0, 1), dim_count(model, :country)))

        # ISM

        if ism_used
            update_param!(inst, :ISMModel_uniforms, rand(Uniform(0, 1), (dim_count(model, :time), dim_count(model, :monsoonsteps))))
        end

        # OMH

        if omh_used
            update_param!(inst, :OMH_uniforms, rand(Uniform(0, 1), dim_count(model, :time)))
        end

        # AMOC

        if amoc_used
            update_param!(inst, :AMOC_uniforms, rand(Uniform(0, 1), dim_count(model, :time)))
        end

        # Amazon

        if amazon_calib != "none"
            update_param!(inst, :AmazonDieback_uniforms, rand(Uniform(0, 1), dim_count(model, :time)))
        end

        # WAIS

        if wais_calib == "Distribution"
            update_param!(inst, :WAISmodel_uniforms, rand(Uniform(0, 1), dim_count(model, :time)))
        end

        run(inst)

        mcres = Dict{Symbol, Any}()
        mcres[:TemperatureModel_T_AT] = copy(inst[:TemperatureModel, :T_AT])
        mcres[:SLRModel_SLR] = copy(inst[:SLRModel, :SLR])
        mcres[:Utility_world_disc_utility] = inst[:Utility, :world_disc_utility]

        if save_rvs
            for jj in 2:ncol(draws)
                mcres[Symbol(names(draws)[jj])] = draws[ii, jj]
            end
        end

        push!(results, mcres)

        next!(progress)
    end

    results
end

function simdataframe(model::Model, results::Vector{Dict{Symbol, Any}}, comp::Symbol, name::Symbol)
    key = Symbol("$(comp)_$(name)")
    if results[1][key] isa Number
        df = DataFrame(trialnum=1:length(results))
        df[!, name] = [results[ii][key] for ii in 1:length(results)]
    else
        dfbase = getdataframe(model, comp, name)
        df = nothing
        for ii in 1:length(results)
            mcdf = dfbase[!, :]
            mcdf[!, name] = results[ii][key]
            mcdf.trialnum .= ii
            if df == nothing
                df = mcdf
            else
                df = vcat(df, mcdf)
            end
        end
    end
    df
end

# function myupdate_param!(inst, unique_name::Symbol, value)
#     if has_parameter(inst.md, unique_name)
#         update_param!(inst, unique_name, value)
#     else
#         parts = split(String(unique_name), "_")
#         set_param!(inst, Symbol(parts[1]), Symbol(join(parts[2:end], "_")), unique_name, value)
#     end
# end
