@defcomp RCP begin
    # Variables
    co2_rcp = Variable(index=[time], unit="GtCO2")
    ch4_conc_rcp = Variable(index=[time], unit="ppb")
    ch4_rcp = Variable(index=[time], unit="MtCH4")
    n2o_conc_rcp = Variable(index=[time], unit="ppb")
    F_EX = Variable(index=[time], unit="W/m^2")
    SO_2 = Variable(index=[time], unit="Mt")
    SO_2_Asia = Variable(index=[time], unit="Mt")

    # ch4_conc_total = Variable(index=[time], unit="ppb")

    # Parameters
    scenario = Parameter{String}()

    # ch4_conc_preind = Parameter(unit="ppb", default=722)
    # ch4_decay = Parameter(default=0.081)
    # ch4_ppb2mt = Parameter(default=2.78)

    # ch4_pcf = Parameter(index=[time], unit="MtCH4")
    # ch4_omh = Parameter(index=[time], unit="MtCH4")
    # ch4_extra = Parameter(index=[time], unit="MtCO2")

    function init(pp, vv, dd)
        scenvals_rma = CSV.read("../data/Scenarios.csv", DataFrame, limit=1, header=0, missingstring="NA")
        scenvals_rcp = CSV.read("../data/RCPs.csv", DataFrame, limit=1, header=0, missingstring="NA")

        scenarios_rma = [scenvals_rma[1, col] for col in names(scenvals_rma)]
        scenarios_rcp = [scenvals_rcp[1, col] for col in names(scenvals_rcp)]

        if pp.scenario ∈ scenarios_rma[2:end]
            rcps = CSV.read("../data/Scenarios.csv", DataFrame, header=2, select=scenarios_rma .== pp.scenario)
        else
            rcps = CSV.read("../data/RCPs.csv", DataFrame, header=2, select=scenarios_rcp .== pp.scenario)
        end
        colnames = replace.(names(rcps), r"_\d+" => "")

        for tt in dd.time
            if gettime(tt) < 2010
                continue
            end
            vv.co2_rcp[TimestepIndex(gettime(tt) - 2009)] = rcps[gettime(tt) - 2009, colnames .== "CO2 (Mt CO2/yr)"][1]
            vv.ch4_rcp[TimestepIndex(gettime(tt) - 2009)] = rcps[gettime(tt) - 2009, colnames .== "CH4 (Mt CH4/yr)"][1]
            vv.ch4_conc_rcp[TimestepIndex(gettime(tt) - 2009)] = rcps[gettime(tt) - 2009, colnames .== "Concentration|CH4 (ppm)"][1]
	    vv.n2o_conc_rcp[TimestepIndex(gettime(tt) - 2009)] = rcps[gettime(tt) - 2009, colnames .== "Concentration|N2O (ppm)"][1]
            vv.F_EX[TimestepIndex(gettime(tt) - 2009)] = rcps[gettime(tt) - 2009, colnames .== "Forcing|other agents (W/m2)"][1]
            vv.SO_2[TimestepIndex(gettime(tt) - 2009)] = rcps[gettime(tt) - 2009, colnames .== "Emissions|Sulfur global (kt SO2/yr)"][1]
            if pp.scenario ∈ scenarios_rma[2:end]
                vv.SO_2_Asia[TimestepIndex(gettime(tt) - 2009)] = rcps[gettime(tt) - 2009, colnames .== "Emissions|Sulfur global (kt SO2/yr)"][1]
            end
        end
    end

    function run_timestep(pp, vv, dd, tt)
        # if is_first(tt)
        #     vv.ch4_conc_total[tt] = pp.ch4_conc_rcp[tt]
        # else
        #     vv.co2_rcp[tt] = pp.co2_rcp[tt] + pp.co2_pcf[tt] + pp.co2_amazon[tt] + pp.co2_extra[tt]
        #     vv.ch4_conc_total[tt] = pp.ch4_conc_preind + (vv.ch4_conc_total[tt-1] - pp.ch4_conc_preind) * (1 - pp.ch4_decay) + (pp.ch4_rcp[tt-1] + pp.ch4_pcf[tt] + pp.ch4_omh[tt])/pp.ch4_ppb2mt
        # end
    end
end

function addRCP(model, scenario; before=nothing, after=nothing)
    scenvals_rma = CSV.read("../data/Scenarios.csv", DataFrame, limit=1, header=0, missingstring="NA")
    scenvals_rcp = CSV.read("../data/RCPs.csv", DataFrame, limit=1, header=0, missingstring="NA")

    scenarios_rma = [scenvals_rma[1, col] for col in names(scenvals_rma)]
    scenarios_rcp = [scenvals_rcp[1, col] for col in names(scenvals_rcp)]
    if (scenario ∉ scenarios_rma[2:end]) && (scenario ∉ scenarios_rcp[2:end])
        throw(ArgumentError("Unknown scenario"))
    end

    emits = add_comp!(model, RCP, first=2010, before=before, after=after)
    emits[:scenario] = scenario

    emits
end
