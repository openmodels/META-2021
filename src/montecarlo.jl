using Mimi
using Distributions
using CSV
using DataFrames
using LinearAlgebra
using ProgressMeter

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
        mcres[:TemperatureModel_T_AT] = inst[:TemperatureModel, :T_AT]
        mcres[:SLRModel_SLR] = inst[:SLRModel, :SLR]
        mcres[:Utility_world_disc_utility] = inst[:Utility, :world_disc_utility]
        mcres[:CO2Model_a0] = inst[:CO2Model, :a0]
        mcres[:CO2Model_a1] = inst[:CO2Model, :a1]
        mcres[:CO2Model_a3] = inst[:CO2Model, :a3]
        mcres[:CO2Model_rho1] = inst[:CO2Model, :rho1]
        mcres[:CO2Model_rho2] = inst[:CO2Model, :rho2]
        mcres[:CO2Model_rho3] = inst[:CO2Model, :rho3]
        mcres[:PostTemperature_r_0] = inst[:PostTemperature, :r_0]
        mcres[:PostTemperature_r_C] = inst[:PostTemperature, :r_C]
        mcres[:PostTemperature_r_T] = inst[:PostTemperature, :r_T]
        mcres[:TemperatureModel_xi_1] = inst[:TemperatureModel, :xi_1]
        mcres[:TemperatureModel_F_2xCO2] = inst[:TemperatureModel, :F_2xCO2]
        mcres[:TemperatureModel_fair_ECS] = inst[:TemperatureModel, :fair_ECS]
        mcres[:TemperatureModel_fair_gamma] = inst[:TemperatureModel, :fair_gamma]
        mcres[:TemperatureModel_fair_C_0] = inst[:TemperatureModel, :fair_C_0]
        mcres[:Forcing_F_2xCO2] = inst[:Forcing, :F_2xCO2]
        mcres[:CH4Model_ch4_alpha] = inst[:CH4Model, :ch4_alpha]
        mcres[:Consumption_damagepersist] = inst[:Consumption, :damagepersist]
        mcres[:Utility_EMUC] = inst[:Utility, :EMUC]
        mcres[:Utility_PRTP] = inst[:Utility, :PRTP]

        push!(results, mcres)

        next!(progress)
    end

    results
end

function getsim(model::Model, pcf_calib::String, amazon_calib::String, gis_calib::String, wais_calib::String, saf_calib::String, ism_used::Bool, omh_used::Bool, amoc_used::Bool, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool; save_rvs::Bool=true)
    sim = @defsim begin

        # CO2Model

        rv(a0_rv) = make_lognormal(0.021089, 0.13636)
        rv(a1) = Logistic(0.224, 0.0000339623)
        rv(a3) = Levy(0.276, 0.00000039629)

        rv(rho1_rv) = make_lognormal(0.0011748, 0.012058)
        rv(rho2_rv) = Beta(1, 1.6911) # Translate to scaled Kumaraswamy in Shifter
        rv(rho3_rv) = make_lognormal(0.017462, 0.29837)

        Shifter.a0_rv = a0_rv
        Shifter.rho1_rv = rho1_rv
        Shifter.rho2_rv = rho2_rv
        Shifter.rho3_rv = rho3_rv

        CO2Model_a1 = a1 # CO2Model.a1
        CO2Model_a3 = a3 # CO2Model.a3

        # PostTemperature

        PostTemperature_r_0 = Normal(32.4, 0.13/1.65*32.4) # PostTemperature.r_0
        PostTemperature_r_C = Normal(0.019, 0.13/1.65*0.019) # PostTemperature.r_C
        PostTemperature_r_T = Normal(4.165, 0.13/1.65*4.165) # PostTemperature.r_T

        # Temperature

        rv(xi_1) = Pareto(5.907, 0.11628)
        rv(F_2xCO2) = Normal(3.45938, 0.43674)
        rv(fair_ECS) = Normal(3.25312, 0.80031)
        rv(fair_gamma) = TriangularDist(0.5, 1.23723, 0.5)
        rv(fair_C_0) = Pareto(1.7062, 53)

        TemperatureModel_xi_1 = xi_1 # TemperatureModel.xi_1
        TemperatureModel_F_2xCO2 = F_2xCO2 # TemperatureModel.F_2xCO2
        TemperatureModel_fair_ECS = fair_ECS # TemperatureModel.fair_ECS
        TemperatureModel_fair_gamma = fair_gamma # TemperatureModel.fair_gamma
        TemperatureModel_fair_C_0 = fair_C_0 # TemperatureModel.fair_C_0

        # Forcing

        Forcing_F_2xCO2 = F_2xCO2 # Forcing.F_2xCO2

        # CH4

        CH4Model_ch4_alpha = TriangularDist(0.0319967, 0.0400033, 0.036) # CH4Model.ch4_alpha

        # SAF

        SAFModel_ECS = fair_ECS # SAFModel.ECS

        # Damages

        Consumption_seeds[:] = DiscreteUniform(1, typemax(Int64)) # Consumption.seeds
        Consumption_slruniforms[:] = Uniform(0, 1) # Consumption.slruniforms

        ## sampling(LHSData, corrlist=all_corrlist)

        save(TemperatureModel.T_AT, SLRModel.SLR,
             Utility.world_disc_utility, CO2Model.a0, CO2Model.a1,
             CO2Model.a3, CO2Model.rho1, CO2Model.rho2, CO2Model.rho3,
             PostTemperature.r_0, PostTemperature.r_C,
             PostTemperature.r_T, TemperatureModel.xi_1,
             TemperatureModel.F_2xCO2, TemperatureModel.fair_ECS,
             TemperatureModel.fair_gamma, TemperatureModel.fair_C_0,
             Forcing.F_2xCO2, CH4Model.ch4_alpha)
    end

    # Interactions

    dists = allinteractrates((symbol, ratemu, ratese) -> Normal(ratemu, ratese))

    add_RV!(sim, :Interactions_gis2amoc, dists["gis"]["amoc"])
    add_RV!(sim, :Interactions_wais2amoc, dists["wais"]["amoc"])
    add_RV!(sim, :Interactions_amaz2amoc, dists["amaz"]["amoc"])
    add_RV!(sim, :Interactions_nino2amoc, dists["nino"]["amoc"])
    add_RV!(sim, :Interactions_amoc2gis, dists["amoc"]["gis"])
    add_RV!(sim, :Interactions_wais2gis, dists["wais"]["gis"])
    add_RV!(sim, :Interactions_amaz2gis, dists["amaz"]["gis"])
    add_RV!(sim, :Interactions_nino2gis, dists["nino"]["gis"])
    add_RV!(sim, :Interactions_amoc2wais, dists["amoc"]["wais"])
    add_RV!(sim, :Interactions_gis2wais, dists["gis"]["wais"])
    add_RV!(sim, :Interactions_amaz2wais, dists["amaz"]["wais"])
    add_RV!(sim, :Interactions_nino2wais, dists["nino"]["wais"])
    add_RV!(sim, :Interactions_amoc2amaz, dists["amoc"]["amaz"])
    add_RV!(sim, :Interactions_gis2amaz, dists["gis"]["amaz"])
    add_RV!(sim, :Interactions_wais2amaz, dists["wais"]["amaz"])
    add_RV!(sim, :Interactions_nino2amaz, dists["nino"]["amaz"])
    add_RV!(sim, :Interactions_amoc2nino, dists["amoc"]["nino"])
    add_RV!(sim, :Interactions_gis2nino, dists["gis"]["nino"])
    add_RV!(sim, :Interactions_wais2nino, dists["wais"]["nino"])
    add_RV!(sim, :Interactions_amaz2nino, dists["amaz"]["nino"])

    if save_rvs
        add_save!(sim, :Interactions, :gis2amoc)
        add_save!(sim, :Interactions, :wais2amoc)
        add_save!(sim, :Interactions, :amaz2amoc)
        add_save!(sim, :Interactions, :nino2amoc)
        add_save!(sim, :Interactions, :amoc2gis)
        add_save!(sim, :Interactions, :wais2gis)
        add_save!(sim, :Interactions, :amaz2gis)
        add_save!(sim, :Interactions, :nino2gis)
        add_save!(sim, :Interactions, :amoc2wais)
        add_save!(sim, :Interactions, :gis2wais)
        add_save!(sim, :Interactions, :amaz2wais)
        add_save!(sim, :Interactions, :nino2wais)
        add_save!(sim, :Interactions, :amoc2amaz)
        add_save!(sim, :Interactions, :gis2amaz)
        add_save!(sim, :Interactions, :wais2amaz)
        add_save!(sim, :Interactions, :nino2amaz)
        add_save!(sim, :Interactions, :amoc2nino)
        add_save!(sim, :Interactions, :gis2nino)
        add_save!(sim, :Interactions, :wais2nino)
        add_save!(sim, :Interactions, :amaz2nino)
    end

    # ISM

    if ism_used
        # ISMModel_uniforms[:, :] = Uniform(0, 1) # ISMModel.uniforms
        add_RV!(sim, :ISMModel_uniforms, ReshapedDistribution([dim_count(model, :time), dim_count(model, :monsoonsteps)], Uniform(0, 1)))
    end

    # OMH

    if omh_used
        # OMH_uniforms[:] = Uniform(0, 1) # OMH.uniforms
        add_RV!(sim, :OMH_uniforms, ReshapedDistribution([dim_count(model, :time)], Uniform(0, 1)))
    end

    # AMOC

    if amoc_used
        # AMOC_uniforms[:] = Uniform(0, 1) # AMOC.uniforms
        add_RV!(sim, :AMOC_uniforms, ReshapedDistribution([dim_count(model, :time)], Uniform(0, 1)))
    end

    # Permafrost

    if pcf_calib == "Kessler probabilistic"
        # PCFModel.propCH4 = Normal(0.023, 0.006)
        # PCFModel.beta_PF = Normal(0.172, 0.0261)
        # PCFModel.C_PF = Normal(1035, 76.53)
        # PCFModel.propPassive = Normal(0.4, 0.055)
        # PCFModel.tau = Normal(70, 30)
        add_RV!(sim, :PCFModel_propCH4, Normal(0.023, 0.006))
        add_RV!(sim, :PCFModel_beta_PF, Normal(0.172, 0.0261))
        add_RV!(sim, :PCFModel_C_PF, Normal(1035, 76.53))
        add_RV!(sim, :PCFModel_propPassive, Normal(0.4, 0.055))
        add_RV!(sim, :PCFModel_tau, Normal(70, 30))

        if save_rvs
            add_save!(sim, :PCFModel, :propCH4)
            add_save!(sim, :PCFModel, :beta_PF)
            add_save!(sim, :PCFModel, :C_PF)
            add_save!(sim, :PCFModel, :propPassive)
            add_save!(sim, :PCFModel, :tau)
        end
    end

    # Amazon

    if amazon_calib == "Distribution"
        # AmazonDieback.Delta_AMAZ = TriangularDist(10, 250, 50)
        add_RV!(sim, :AmazonDieback_Delta_AMAZ, TriangularDist(10, 250, 50))

        if save_rvs
            add_save!(sim, :AmazonDieback, :Delta_AMAZ)
        end
    end
    if amazon_calib != "none"
        # AmazonDieback_uniforms[:] = Uniform(0, 1) # AmazonDieback.uniforms
        add_RV!(sim, :AmazonDieback_uniforms, ReshapedDistribution([dim_count(model, :time)], Uniform(0, 1)))
    end

    # GIS

    if gis_calib == "Distribution"
        # GISModel.avoldot = Normal(-0.0000106, 0.0000244/0.5/100) # Only works with a meltmult of 1
        add_RV!(sim, :GISModel_avoldot, Normal(-0.0000106, 0.0000244/0.5/100)) # Only works with a meltmult of 1

        if save_rvs
            add_save!(sim, :GISModel, :avoldot)
        end
    end

    # WAIS

    if wais_calib == "Distribution"
        # WAISmodel.waisrate = make_lognormal(3.3 / 1000, 1.65 / 1000)
        add_RV!(sim, :WAISmodel_waisrate, make_lognormal(3.3 / 1000, 1.65 / 1000))
        # WAISmodel_uniforms[:] = Uniform(0, 1) # WAISmodel.uniforms
        add_RV!(sim, :WAISmodel_uniforms, [Uniform(0, 1) for tt in 1:dim_count(model, :time)])

        if save_rvs
            add_save!(sim, :WAISmodel, :waisrate)
        end
    end

    # SAF

    if saf_calib == "Distribution"
        # SAFModel.delta = TriangularDist(-1, 1, 0)
        # SAFModel.FRT = TriangularDist(10, 55, 20)
        add_RV!(sim, :SAFModel_saf_delta, TriangularDist(-1, 1, 0))
        add_RV!(sim, :SAFModel_FRT, TriangularDist(10, 55, 20))

        if save_rvs
            add_save!(sim, :SAFModel, :saf_delta)
            add_save!(sim, :SAFModel, :FRT)
        end
    end

    # Damages

    if persist_dist
        add_RV!(sim, :Consumption_damagepersist, Uniform(0, 1))

        if save_rvs
            add_save!(sim, :Consumption, :damagepersist)
        end
    end

    # Utility

    if emuc_dist
        add_RV!(sim, :Consumption_EMUC, TriangularDist(0.5, 2, 1.5))
        if save_rvs
            add_save!(sim, :Consumption, :EMUC)
        end
    end
    if prtp_dist
        add_RV!(sim, :Consumption_PRTP, TriangularDist(0.001, 0.02, 0.01))
        if save_rvs
            add_save!(sim, :Consumption, :PRTP)
        end
    end

    sim
end
