@defcomp OMH begin
    # Variables
    p_OMH = Variable(index=[time])
    I_OMH = Variable{Bool}(index=[time])
    CH4_OMH = Variable(index=[time], unit="MtCH4")
    cum_CH4_OMH = Variable(index=[time], unit="MtCH4")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC")

    uniforms = Parameter(index=[time])
    b_OMH = Parameter(unit="1/degC")
    max_CH4_OMH = Parameter(unit="MtCH4")
    Delta_OMH = Parameter(unit="year")

    function run_timestep(pp, vv, dd, tt)
        # Calculate probability of triggering
        vv.p_OMH[tt] = 1 - exp(-pp.b_OMH * pp.T_AT[tt])

        ## Determine if OMH is triggered
        if is_first(tt)
            vv.I_OMH[tt] = pp.uniforms[tt] < vv.p_OMH[tt]

            # Calculate CH4 release
            vv.CH4_OMH[tt] = vv.I_OMH[tt] ? pp.max_CH4_OMH / pp.Delta_OMH : 0
            vv.cum_CH4_OMH[tt] = vv.CH4_OMH[tt]
        else
            vv.I_OMH[tt] = vv.I_OMH[tt-1] || (pp.uniforms[tt] < vv.p_OMH[tt])

            # Calculate CH4 release
            vv.CH4_OMH[tt] = (vv.I_OMH[tt] && vv.cum_CH4_OMH[tt-1] < pp.max_CH4_OMH) ?
                pp.max_CH4_OMH / pp.Delta_OMH : 0
            vv.cum_CH4_OMH[tt] = vv.cum_CH4_OMH[tt-1] + vv.CH4_OMH[tt]
        end
    end
end

function addOMH(model, calibration; before=nothing, after=nothing)
    params = CSV.read("../data/OMH.csv", DataFrame)
    if calibration ∉ params.Calibration
        throw(ArgumentError("Unknown OMH calibration"))
    end

    omh = add_comp!(model, OMH, first=2010, before=before, after=after)
    omh[:b_OMH] = params[params.Calibration .== calibration, "b_OMH"][1]
    omh[:max_CH4_OMH] = params[params.Calibration .== calibration, "CH4_OMH max (MtCH4)"][1]
    omh[:Delta_OMH] = params[params.Calibration .== calibration, "Delta_OMH"][1]

    omh
end
