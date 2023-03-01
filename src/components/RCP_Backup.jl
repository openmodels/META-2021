@defcomp RCP begin
    # Variables
    co2_rcp = Variable(index=[time], unit="GtCO2")
    ch4_conc_rcp = Variable(index=[time], unit="ppb")
    ch4_rcp = Variable(index=[time], unit="MtCH4")

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
        scenarios = vec(convert(Array, CSV.read("../data/RCPs.csv", DataFrame, limit=1, header=0, missingstring="NA")))
        rcps = CSV.read("../data/RCPs.csv", DataFrame, header=2, select=scenarios .== scenario)
        colnames = replace.(names(rcps), r"_\d+" => "")

        for tt in 1:dd.time
            vv.co2_rcp[TimeIndex(tt)] = rcps[!, colnames .== "CO2 (Mt CO2/yr)"][!, 1]
            vv.ch4_rcp[TimeIndex(tt)] = rcp[!, colnames .== "CH4 (Mt CH4/yr)"][!, 1]
            vv.ch4_conc_rcp[TimeIndex(tt)] = rcp[!, colnames .== "Concentration|CH4"][!, 1]
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

function addRCP(model, scenario)
    scenarios = vec(convert(Array, CSV.read("../data/RCPs.csv", DataFrame, limit=1, header=0, missingstring="NA")))
    if scenario ∉ scenarios[2:end]
        throw(ArgumentError("Unknown RCP scenario"))
    end

    emits = add_comp!(model, Emissions)
    emits[:scenario] = scenario

    emits
end
