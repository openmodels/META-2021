using Mimi
using Distributions
using CSV
using DataFrames

@defcomp Shifter begin
    a0 = Variable()
    rho1 = Variable()
    rho2 = Variable()
    rho3 = Variable()

    a0_rv = Parameter()
    rho1_rv = Parameter()
    rho2_rv = Parameter()
    rho3_rv = Parameter()

    function init(pp, vv, dd)
        vv.a0 = pp.a0_rv + 0.213777
        vv.rho1 = pp.rho1_rv + 0.0024261
        vv.rho2 = (0.048002 - 0.01674) * pp.rho2_rv^(1 / 1.013) + 0.01674
        vv.rho3 = pp.rho3_rv + 0.231322
    end
end

# Prepare correlations for CO2 model
## Note: Correlations with rho2 are applied to Beta, not Kumaraswamy
all_corrlist = []
function add_corrmat(corrs::DataFrame, skips::Vector{String})
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

function getsim(permafrost_calib, amazon_calib, gis_calib, wais_calib, saf_calib)
    sim = @defsim begin

        ## CO2Model

        rv(a0_rv) = LogNormal(0.021089, 0.13636)
        rv(a1) = Logistic(0.224, 0.0000339623)
        rv(a3) = Levy(0.276, 0.00000039629)

        rv(rho1_rv) = LogNormal(0.0011748, 0.012058)
        rv(rho2_rv) = Beta(1, 1.6911) # Translate to scaled Kumaraswamy in Shifter
        rv(rho3_rv) = LogNormal(0.017462, 0.29837)

        # PostTemperature

        PostTemperature.r_0 = Normal(32.4, 0.13/1.65*32.4)
        PostTemperature.r_C = Normal(0.019, 0.13/1.65*0.019)
        PostTemperature.r_T = Normal(4.165, 0.13/1.65*4.165)

        # Temperature

        rv(xi_1) = Pareto(0.11628, 5.907)
        rv(F_2xCO2) = Normal(3.45938, 0.43674)
        rv(ECS) = Normal(3.25312, 0.80031)
        rv(gamma) = TriangularDist(0.5, 0.5, 1.23723)
        rv(C_0) = Pareto(53, 1.7062)

        # CH4

        CH4Model.alpha = TriangularDist(0.0319967, 0.036, 0.0400033)

        # Permafrost

        if permafrost_calib == "Kessler probabilistic"
            PCFModel.propCH4 = Normal(0.023, 0.006)
            PCFModel.beta_PF = Normal(0.172, 0.0261)
            PCFModel.C_PF = Normal(1035, 76.53)
            PCFModel.propPassive = Normal(0.4, 0.055)
            PCFModel.tau = Normal(70, 30)
        end

        # Amazon

        if amazon_calib == "Distribution"
            AmazonDieback.Delta_AMAZ = TriangularDist(10, 50, 250)
        end

        # GIS

        if gis_calib == "Distribution"
            GISModel.avoldot = Normal(-0.0000106, 0.0000244/0.5/100) # Only works with a meltmult of 1
        end

        # WAIS

        if wais_calib == "Distribution"
            WAISmodel.waisrate = LogNormal(3.3 / log(1000), 1.65 / log(1000))
        end

        # SAF

        if saf_calib == "Distribution"
            SAFModel.delta = TriangularDist(-1, 0, 1)
            SAFModel.FRT = TriangularDist(10, 20, 55)
        end

        # Interactions

        dists = allinteractrates((symbol, ratemu, ratese) -> Normal(ratemu, ratese))
        Interactions.gis2amoc = dists["gis"]["amoc"]
        Interactions.wais2amoc = dists["wais"]["amoc"]
        Interactions.amaz2amoc = dists["amaz"]["amoc"]
        Interactions.nino2amoc = dists["nino"]["amoc"]
        Interactions.amoc2gis = dists["amoc"]["gis"]
        Interactions.wais2gis = dists["wais"]["gis"]
        Interactions.amaz2gis = dists["amaz"]["gis"]
        Interactions.nino2gis = dists["nino"]["gis"]
        Interactions.amoc2wais = dists["amoc"]["wais"]
        Interactions.gis2wais = dists["gis"]["wais"]
        Interactions.amaz2wais = dists["amaz"]["wais"]
        Interactions.nino2wais = dists["nino"]["wais"]
        Interactions.amoc2amaz = dists["amoc"]["amaz"]
        Interactions.gis2amaz = dists["gis"]["amaz"]
        Interactions.wais2amaz = dists["wais"]["amaz"]
        Interactions.nino2amaz = dists["nino"]["amaz"]
        Interactions.amoc2nino = dists["amoc"]["nino"]
        Interactions.gis2nino = dists["gis"]["nino"]
        Interactions.wais2nino = dists["wais"]["nino"]
        Interactions.amaz2nino = dists["amaz"]["nino"]

        # Damages

        Consumption.seeds = DiscreteUniform(1, typemax(Int64))
        Consumption.damagepersist = Uniform(0, 1)
        Consumption.slruniforms = Uniform(0, 1)

        # Utility
        Utility.EMUC = TriangularDist(0.5, 2, 1.5)
        Utility.PRTP = TriangularDist(0.001, 0.02, 0.01)

        sampling(LHSData, corrlist=all_corrlist)

        # # Indicate which variables to save for each model run.
        # # The syntax is: component_name.variable_name
        # save(grosseconomy.K, grosseconomy.YGROSS,
        #      emissions.E, emissions.E_Global)
    end

    sim
end
