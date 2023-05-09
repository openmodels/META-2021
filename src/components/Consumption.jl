include("../lib/gdppc.jl")
include("../lib/damages.jl")
include("../lib/saverate.jl")

@defcomp Consumption begin
    # Variables
    gdppc_region = Variable(index=[time, region], unit="2010 USD PPP")
    gdppc_ratio_region = Variable(index=[time, region])
    gdppc_growth_region = Variable(index=[time, region])

    gdppc_growth = Variable(index=[time, country])
    gdppc = Variable(index=[time, country], unit="2010 USD PPP")
    conspc_preadj = Variable(index=[time, country], unit="2010 USD PPP") # previous year's
    conspc = Variable(index=[time, country], unit="2010 USD PPP")
    baseline_consumption_percap_percountry = Variable(index = [time, country], unit = "2010 USD PPP") # Counterfactual consumption per cap per country from SSPs

    # Parameters
    ssp = Parameter{String}()

    # Based on SSP
    convergerate_gdppc = Parameter()
    decayrate_gdppc = Parameter()
    popweights_region = Parameter(index=[region])

    # Based on historical data
    gdppc_2009 = Parameter(index=[country], unit="2010 USD PPP")
    saverate = Parameter(index=[country])

    # Based on damage specification
    beta1 = Parameter(unit="1/degC")
    beta2 = Parameter(unit="1/degC^2")

    slrdamageconfig = Parameter{String}()
    slruniforms = Parameter(index=[country]) # only used for MC mode
    slrcoeff = Parameter(index=[country], unit="1/m")

    # Damage configuration
    damagepersist = Parameter(default=0.5)
    min_conspc = Parameter(unit="2010 USD PPP", default=1)
    extradamage = Parameter(index=[time, country])

    # Inputs from other components
    T_country = Parameter(index=[time, country], unit="degC")
    SLR = Parameter(index=[time], unit="m")

    function init(pp, vv, dd)
        isos = dim_keys(model, :country)

        if pp.slrdamageconfig == "distribution"
            pp.slrcoeff = [getslrcoeff_distribution(isos[cc], slrdamage, pp.slruniforms[cc]) for cc in 1:length(isos)]
        end
    end

    function run_timestep(pp, vv, dd, tt)
        isos = dim_keys(model, :country)

        if is_first(tt)
            for rr in dd.region
                vv.gdppc_region[tt, rr] = getgdppc_ssp(dim_keys(model, :region)[rr], pp.ssp, 2010)
                if ismissing(vv.gdppc_region[tt, rr])
                    vv.gdppc_region[tt, rr] = 0
                end
                vv.gdppc_ratio_region[tt, rr] = 1
                vv.gdppc_growth_region[tt, rr] = 0 # ignored
            end

            for cc in dd.country
                gdppc = getgdppc(isos[cc], 2010)
                if ismissing(gdppc)
                    if isos[cc] == "DJI"
                        vv.gdppc[tt, cc] = getgdppc(isos[cc], 2011)
                        vv.gdppc_growth[tt, cc] = vv.gdppc[tt, cc] / pp.gdppc_2009[cc] - 1
                    else
                        vv.gdppc[tt, cc] = vv.gdppc_growth[tt, cc] = 0
                    end
                else
                    vv.gdppc[tt, cc] = gdppc
                    vv.gdppc_growth[tt, cc] = vv.gdppc[tt, cc] / pp.gdppc_2009[cc] - 1
                end

                vv.conspc_preadj[tt, cc] = (1-pp.saverate[cc])*pp.gdppc_2009[cc]
            end
        else
            for rr in dd.region
                vv.gdppc_region[tt, rr] = getgdppc_ssp(dim_keys(model, :region)[rr], pp.ssp, gettime(tt))
                vv.gdppc_ratio_region[tt, rr] = vv.gdppc_region[tt, rr] / vv.gdppc_region[TimestepIndex(1), rr]
                if gettime(tt) <= 2100
                    vv.gdppc_growth_region[tt, rr] = vv.gdppc_ratio_region[tt, rr] / vv.gdppc_ratio_region[tt-1, rr] - 1
                else
                    vv.gdppc_growth_region[tt, rr] = (1-pp.convergerate_gdppc-pp.decayrate_gdppc)*vv.gdppc_growth_region[tt-1, rr]+pp.decayrate_gdppc*sum(vv.gdppc_growth_region[tt-1, :] .* pp.popweights_region)
                end
            end

            for cc in dd.country
                region = getregion(isos[cc])
                if ismissing(region)
                    growth = 0
                else
                    rr = findfirst(dim_keys(model, :region) .== region)
                    growth = vv.gdppc_growth_region[tt, rr]
                end

                vv.gdppc_growth[tt, cc] = growth
                vv.gdppc[tt, cc] = vv.gdppc[tt-1, cc] * (1 + growth)

                if vv.conspc[tt-1, cc] == 0 || vv.gdppc[tt, cc] == 0
                    vv.conspc_preadj[tt, cc] = 0
                else
                    vv.conspc_preadj[tt, cc] = pp.damagepersist*(1-pp.saverate[cc])*vv.gdppc[tt-1, cc]+(1-pp.damagepersist)*vv.conspc[tt-1, cc]
                end
            end
        end

        for cc in dd.country
            vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+(vv.gdppc_growth[tt, cc]+pp.beta1*(pp.T_country[tt, cc]-pp.T_country[TimestepIndex(1), cc])+pp.beta2*(pp.T_country[tt, cc]^2-pp.T_country[TimestepIndex(1), cc]^2)))*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])

            # Compute baseline consumption per capita without damages
            vv.baseline_consumption_percap_percountry[tt,cc] = (1-pp.saverate[cc])*vv.gdppc[tt, cc]

            if vv.conspc[tt, cc] <= pp.min_conspc && vv.conspc[TimestepIndex(1), cc] != 0
                vv.conspc[tt, cc] = pp.min_conspc
            end
        end
    end
end

function addConsumption(model, tdamage, slrdamage, ssp)
    if tdamage ∉ ["none", "distribution", "pointestimate", "low", "high"]
        throw(ArgumentError("Unknown Consumption tdamage"))
    end
    if slrdamage ∉ ["none", "distribution", "mode", "low", "high"]
        throw(ArgumentError("Unknown Consumption slrdamage"))
    end

    cons = add_comp!(model, Consumption, first=2010)

    cons[:ssp] = ssp

    sspextend = CSV.read("../data/sspextend-gdppc.csv", DataFrame)
    cons[:convergerate_gdppc] = sspextend[sspextend.SSP .== ssp, "Convergence rate"][1]
    cons[:decayrate_gdppc] = sspextend[sspextend.SSP .== ssp, "Decay rate"][1]

    cons[:popweights_region] = [getpopweight_ssp(region, ssp) for region in dim_keys(model, :region)]

    isos = dim_keys(model, :country)

    if tdamage == "none"
        cons[:beta1] = cons[:beta2] = 0.
    else
        cons[:beta1], cons[:beta2] = getbhmbetas(tdamage)
    end

    if slrdamage == "none"
        cons[:slrcoeff] = zeros(length(isos))
        cons[:slruniforms] = zeros(length(isos))
    elseif slrdamage != "distribution"
        cons[:slrcoeff] = [getslrcoeff(iso, slrdamage) for iso in isos]
        cons[:slruniforms] = zeros(length(isos))
    else
        cons[:slrcoeff] = zeros(length(isos))  # to be filled by init
        cons[:slruniforms] = zeros(length(isos))  # to be filled by monte carlo
    end

    cons[:slrdamageconfig] = slrdamage

    cons[:saverate] = [getsaverate(iso) for iso in isos]
    gdppc_2009 = [getgdppc(iso, 2009) for iso in isos]
    gdppc_2009[ismissing.(gdppc_2009)] .= 0
    gdppc_2009[isos .== "DJI"] .= 2700
    cons[:gdppc_2009] = gdppc_2009
    cons[:extradamage] = zeros(dim_count(model, :time), dim_count(model, :country))

    cons
end

