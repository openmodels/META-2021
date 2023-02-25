using Mimi
using Distributions
using CSV
using DataFrames
using LinearAlgebra
using ProgressMeter
include("../src/lib/presets.jl")

import Mimi.add_RV!, Mimi.add_save!, Mimi.ReshapedDistribution, Mimi.lhs, Mimi.RandomVariable

# Prepare correlations for CO2 model
## Note: Correlations with rho2 are applied to Beta, not Kumaraswamy
function load_corrmat(filepath::String, skips::Vector{String}=String[])
    corrs = CSV.read(filepath, DataFrame)
    corrmatrix = Matrix(1.0I, nrow(corrs), nrow(corrs))

    for ii in 2:nrow(corrs)
        if corrs[ii, 1] in skips
            continue
        end
        for jj in 2:ii
            if names(corrs)[jj] in skips
                continue
            end
            corrmatrix[ii, jj-1] = corrmatrix[jj-1, ii] = corrs[ii, jj]
        end
    end

    if skips != []
        iis = []
        for ii in 1:nrow(corrs)
            if corrs[ii, 1] ∉ skips
                push!(iis, ii)
            end
        end
        corrmatrix = corrmatrix[iis, iis]
    end

    corrmatrix
end

function make_lognormal(riskmu, risksd)
    mu = log(riskmu^2 / sqrt(risksd^2 + riskmu^2))
    sd = sqrt(log(1 + (risksd /  riskmu)^2))
    LogNormal(mu, sd)
end

# abstract type SimulationDefAbstract end

# mutable struct FunctionalSimulationDef <: SimulationDefAbstract
#     mcupdate::Function
# end


function getsim_base(trials::Int64, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool)
    draws = DataFrame(mc=1:trials)
    # CO2Model

    corrmatrix = load_corrmat("../data/arhocorrs.csv", ["a2"])
    rvlist = Vector{RandomVariable}([RandomVariable(:a0_rv, make_lognormal(0.021089, 0.13636)),
                                     RandomVariable(:a1, Logistic(0.224, 0.0000339623)),
                                     RandomVariable(:a3, Levy(0.276, 0.00000039629)),
                                     RandomVariable(:rho1_rv, make_lognormal(0.0011748, 0.012058)),
                                     RandomVariable(:rho2_rv, Beta(1, 1.6911)), # Translate to scaled Kumaraswamy
                                     RandomVariable(:rho3_rv, make_lognormal(0.017462, 0.29837))])
    df = lhs(rvlist, trials; corrmatrix=corrmatrix)

    draws.CO2Model_a0 = df.a0_rv .+ 0.213777
    draws.CO2Model_a1 = df.a1
    draws.CO2Model_a3 = df.a3

    draws.CO2Model_rho1 = df.rho1_rv .+ 0.0024261
    draws.CO2Model_rho2 = (0.048002 - 0.01674) * df.rho2_rv.^(1 / 1.013) .+ 0.01674
    draws.CO2Model_rho3 = df.rho3_rv .+ 0.231322

    # PostTemperature

    draws.PostTemperature_r_0 = rand(Normal(32.4, 0.13/1.65*32.4), trials) # PostTemperature.r_0
    draws.PostTemperature_r_C = rand(Normal(0.019, 0.13/1.65*0.019), trials) # PostTemperature.r_C
    draws.PostTemperature_r_T = rand(Normal(4.165, 0.13/1.65*4.165), trials) # PostTemperature.r_T

    # Temperature

    corrmatrix = load_corrmat("../data/tempcorrs.csv")
    rvlist = Vector{RandomVariable}([RandomVariable(:xi_1, Pareto(5.907, 0.11628)),
                                     RandomVariable(:F_2xCO2,  Normal(3.45938, 0.43674)),
                                     RandomVariable(:fair_ECS,  Normal(3.25312, 0.80031)),
                                     RandomVariable(:fair_gamma,  TriangularDist(0.5, 1.23723, 0.5)),
                                     RandomVariable(:fair_C_0,  Pareto(1.7062, 53))])
    df = lhs(rvlist, trials; corrmatrix=corrmatrix)

    draws.TemperatureModel_xi_1 = df.xi_1 # TemperatureModel.xi_1
    draws.TemperatureModel_F_2xCO2 = df.F_2xCO2 # TemperatureModel.F_2xCO2
    draws.TemperatureModel_fair_ECS = df.fair_ECS # TemperatureModel.fair_ECS
    draws.TemperatureModel_fair_gamma = df.fair_gamma # TemperatureModel.fair_gamma
    draws.TemperatureModel_fair_C_0 = df.fair_C_0 # TemperatureModel.fair_C_0

    # Forcing

    draws.Forcing_F_2xCO2 = df.F_2xCO2 # Forcing.F_2xCO2

    # CH4

    draws.CH4Model_ch4_alpha = rand(TriangularDist(0.0319967, 0.0400033, 0.036), trials) # CH4Model.ch4_alpha

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

function runsim_base(model::Model, draws::DataFrame; save_rvs::Bool=true)
    inst = Mimi.build(model)

    progress = Progress(nrow(draws), 1)

    results = Dict{Symbol, Any}[]
    for ii in 1:nrow(draws)
        for jj in 2:ncol(draws)
            update_param!(inst, Symbol(names(draws)[jj]), draws[ii, jj])
        end

        # Damages

        update_param!(inst, :Consumption_seeds, rand(DiscreteUniform(1, typemax(Int64)), dim_count(model, :country)))
        update_param!(inst, :Consumption_slruniforms, rand(Uniform(0, 1), dim_count(model, :country)))

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

function getsim(trials::Int64, pcf_calib::String, amazon_calib::String, gis_calib::String, wais_calib::String, saf_calib::String, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool)
    draws = getsim_base(trials, persist_dist, emuc_dist, prtp_dist)

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

    inst = Mimi.build(model)

    progress = Progress(nrow(draws), 1)

    results = Dict{Symbol, Any}[]
    for ii in 1:nrow(draws)
        for jj in 2:ncol(draws)
            if has_parameter(model.md, Symbol(names(draws)[jj]))
                update_param!(inst, Symbol(names(draws)[jj]), draws[ii, jj])
            end
        end

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
