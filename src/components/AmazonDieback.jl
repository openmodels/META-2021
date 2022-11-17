@defcomp AmazonDieback begin
    # Variables
    p_AMAZ = Variable(index=[time])
    I_AMAZ = Variable{Bool}(index=[time])
    CO2_AMAZ = Variable(index=[time], unit="GtCO2")
    cum_CO2_AMAZ = Variable(index=[time], unit="GtCO2")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC")
    probmult = Parameter(index=[time])

    uniforms = Parameter(index=[time])
    b_AMAZ = Parameter(unit="1/degC")
    max_CO2_AMAZ = Parameter(unit="GtCO2")
    Delta_AMAZ = Parameter(unit="year")


    function run_timestep(pp, vv, dd, tt)
        # Calculate probability of triggering
        vv.p_AMAZ[tt] = min((1 - exp(-pp.b_AMAZ * max(0, pp.T_AT[tt]-1)))*pp.probmult[tt], 1)

        ## Determine if Amazon dieback is triggered
        if is_first(tt)
            vv.I_AMAZ[tt] = pp.uniforms[tt] < vv.p_AMAZ[tt]

            # Calculate CO2 release
            vv.CO2_AMAZ[tt] = vv.I_AMAZ[tt] ? pp.max_CO2_AMAZ / pp.Delta_AMAZ : 0
            vv.cum_CO2_AMAZ[tt] = vv.CO2_AMAZ[tt]
        else
            vv.I_AMAZ[tt] = vv.I_AMAZ[tt-1] || (pp.uniforms[tt] < vv.p_AMAZ[tt])

            # Calculate CO2 release
            vv.CO2_AMAZ[tt] = (vv.I_AMAZ[tt] && vv.cum_CO2_AMAZ[tt-1] < pp.max_CO2_AMAZ - 0.01) ?
                pp.max_CO2_AMAZ / pp.Delta_AMAZ : 0
            vv.cum_CO2_AMAZ[tt] = vv.cum_CO2_AMAZ[tt-1] + vv.CO2_AMAZ[tt]
        end
    end
end

function addAmazonDieback(model, calibration; before=nothing, after=nothing)
    params = CSV.read("../data/AMAZparams.csv", DataFrame)
    if calibration ∉ params.Calibration
        throw(ArgumentError("Unknown AMAZ calibration"))
    end

    amazon = add_comp!(model, AmazonDieback, before=before, after=after)
    amazon[:b_AMAZ] = params[params.Calibration .== calibration, "b_AMAZ"][1]
    amazon[:max_CO2_AMAZ] = params[params.Calibration .== calibration, "CO2_AMAZ max (GtCO2)"][1]
    amazon[:Delta_AMAZ] = params[params.Calibration .== calibration, "Delta_AMAZ"][1]

    amazon[:probmult] = ones(dim_count(model, :time))

    amazon
end
