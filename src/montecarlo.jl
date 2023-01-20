using Mimi
using Distributions
using CSV
using DataFrames

import Mimi.add_RV!

@defcomp Shifter begin
    a0 = Variable()
    rho1 = Variable()
    rho2 = Variable()
    rho3 = Variable()

    a0_rv = Parameter(default=0.)
    rho1_rv = Parameter(default=0.)
    rho2_rv = Parameter(default=0.)
    rho3_rv = Parameter(default=0.)

    function init(pp, vv, dd)
        vv.a0 = pp.a0_rv + 0.213777
        vv.rho1 = pp.rho1_rv + 0.0024261
        vv.rho2 = (0.048002 - 0.01674) * pp.rho2_rv^(1 / 1.013) + 0.01674
        vv.rho3 = pp.rho3_rv + 0.231322
    end
end

# Prepare correlations for CO2 model
## Note: Correlations with rho2 are applied to Beta, not Kumaraswamy
all_corrlist = Tuple{Symbol, Symbol, Float64}[]
function add_corrmat(corrs::DataFrame, skips::Vector{String}=String[])
    for ii in 2:nrow(corrs)
        if corrs[ii, 1] in skips
            continue
        end
        for jj in 2:ii
            if names(corrs)[jj] in skips
                continue
            end
            push!(all_corrlist, (Symbol(names(corrs)[jj]), Symbol(corrs[ii, 1]), corrs[ii, jj]))
        end
    end
end

arhocorrs = CSV.read("../data/arhocorrs.csv", DataFrame)
add_corrmat(arhocorrs, ["a2"])

tempcorrs = CSV.read("../data/tempcorrs.csv", DataFrame)
add_corrmat(tempcorrs)

function prepare_montecarlo!(model::Model)
    add_comp!(model, Shifter, before=:CO2Model)

    connect_param!(model, :CO2Model, :a0, :Shifter, :a0)
    connect_param!(model, :CO2Model, :rho1, :Shifter, :rho1)
    connect_param!(model, :CO2Model, :rho2, :Shifter, :rho2)
    connect_param!(model, :CO2Model, :rho3, :Shifter, :rho3)
end

function make_lognormal(riskmu, risksd)
    mu = log(riskmu^2 / sqrt(risksd^2 + riskmu^2))
    sd = sqrt(log(1 + (risksd /  riskmu)^2))
    LogNormal(mu, sd)
end

function getsim(pcf_calib::String, amazon_calib::String, gis_calib::String, wais_calib::String, saf_calib::String, persist_dist::Bool, emuc_dist::Bool, prtp_dist::Bool)
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
        rv(fair_C_0) = Pareto(53, 1.7062)

        TemperatureModel_xi_1 = xi_1 # TemperatureModel.xi_1
        TemperatureModel_F_2xCO2 = F_2xCO2 # TemperatureModel.F_2xCO2
        TemperatureModel_fair_ECS = fair_ECS # TemperatureModel.fair_ECS
        TemperatureModel_fair_gamma = fair_gamma # TemperatureModel.fair_gamma
        TemperatureModel_fair_C_0 = fair_C_0 # TemperatureModel.fair_C_0

        # Forcing

        Forcing_F_2xCO2 = F_2xCO2 # Forcing.F_2xCO2

        # CH4

        CH4Model_ch4_alpha = TriangularDist(0.0319967, 0.0400033, 0.036) # CH4Model.ch4_alpha

        # ISM

        ISMModel_uniforms[:, :] = Uniform(0, 1) # ISMModel.uniforms

        # OMH

        OMH_uniforms[:] = Uniform(0, 1) # OMH.uniforms

        # Amazonn Dieback

        AmazonDieback_uniforms[:] = Uniform(0, 1) # AmazonDieback.uniforms

        # WAIS

        WAISmodel_uniforms[:] = Uniform(0, 1) # WAISmodel.uniforms

        # AMOC

        AMOC_uniforms[:] = Uniform(0, 1) # AMOC.uniforms

        # SAF

        SAFModel_ECS = fair_ECS # SAFModel.ECS

        # Damages

        Consumption_seeds[:] = DiscreteUniform(1, typemax(Int64)) # Consumption.seeds
        Consumption_slruniforms[:] = Uniform(0, 1) # Consumption.slruniforms

        sampling(LHSData, corrlist=all_corrlist)

        save(TemperatureModel.T_AT, SLRModel.SLR, Utility.world_disc_utility)
    end

    # Interactions

    dists = allinteractrates((symbol, ratemu, ratese) -> Normal(ratemu, ratese))

    add_RV!(sim, :gis2amoc, dists["gis"]["amoc"])
    add_RV!(sim, :wais2amoc, dists["wais"]["amoc"])
    add_RV!(sim, :amaz2amoc, dists["amaz"]["amoc"])
    add_RV!(sim, :nino2amoc, dists["nino"]["amoc"])
    add_RV!(sim, :amoc2gis, dists["amoc"]["gis"])
    add_RV!(sim, :wais2gis, dists["wais"]["gis"])
    add_RV!(sim, :amaz2gis, dists["amaz"]["gis"])
    add_RV!(sim, :nino2gis, dists["nino"]["gis"])
    add_RV!(sim, :amoc2wais, dists["amoc"]["wais"])
    add_RV!(sim, :gis2wais, dists["gis"]["wais"])
    add_RV!(sim, :amaz2wais, dists["amaz"]["wais"])
    add_RV!(sim, :nino2wais, dists["nino"]["wais"])
    add_RV!(sim, :amoc2amaz, dists["amoc"]["amaz"])
    add_RV!(sim, :gis2amaz, dists["gis"]["amaz"])
    add_RV!(sim, :wais2amaz, dists["wais"]["amaz"])
    add_RV!(sim, :nino2amaz, dists["nino"]["amaz"])
    add_RV!(sim, :amoc2nino, dists["amoc"]["nino"])
    add_RV!(sim, :gis2nino, dists["gis"]["nino"])
    add_RV!(sim, :wais2nino, dists["wais"]["nino"])
    add_RV!(sim, :amaz2nino, dists["amaz"]["nino"])

    # Permafrost

    if pcf_calib == "Kessler probabilistic"
        # PCFModel.propCH4 = Normal(0.023, 0.006)
        # PCFModel.beta_PF = Normal(0.172, 0.0261)
        # PCFModel.C_PF = Normal(1035, 76.53)
        # PCFModel.propPassive = Normal(0.4, 0.055)
        # PCFModel.tau = Normal(70, 30)
        add_RV!(sim, :propCH4, Normal(0.023, 0.006))
        add_RV!(sim, :beta_PF, Normal(0.172, 0.0261))
        add_RV!(sim, :C_PF, Normal(1035, 76.53))
        add_RV!(sim, :propPassive, Normal(0.4, 0.055))
        add_RV!(sim, :tau, Normal(70, 30))
    end

    # Amazon

    if amazon_calib == "Distribution"
        # AmazonDieback.Delta_AMAZ = TriangularDist(10, 250, 50)
        add_RV!(sim, :Delta_AMAZ, TriangularDist(10, 250, 50))
    end

    # GIS

    if gis_calib == "Distribution"
        # GISModel.avoldot = Normal(-0.0000106, 0.0000244/0.5/100) # Only works with a meltmult of 1
        add_RV!(sim, :avoldot, Normal(-0.0000106, 0.0000244/0.5/100)) # Only works with a meltmult of 1
    end

    # WAIS

    if wais_calib == "Distribution"
        # WAISmodel.waisrate = make_lognormal(3.3 / log(1000), 1.65 / log(1000))
        add_RV!(sim, :waisrate, make_lognormal(3.3 / log(1000), 1.65 / log(1000)))
    end

    # SAF

    if saf_calib == "Distribution"
        # SAFModel.delta = TriangularDist(-1, 1, 0)
        # SAFModel.FRT = TriangularDist(10, 55, 20)
        add_RV!(sim, :delta, TriangularDist(-1, 1, 0))
        add_RV!(sim, :FRT, TriangularDist(10, 55, 20))
    end

    # Damages

    if persist_dist
        add_RV!(sim, :damagepersist, Uniform(0, 1))
    end

    # Utility

    if emuc_dist
        add_RV!(sim, :EMUC, TriangularDist(0.5, 2, 1.5))
    end
    if prtp_dist
        add_RV!(sim, :PRTP, TriangularDist(0.001, 0.02, 0.01))
    end

    sim
end
