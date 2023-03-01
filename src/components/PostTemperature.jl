@defcomp PostTemperature begin
    # Variables
    iIRF_100 = Variable(index=[time], unit="GtCO2")
    alpha = Variable(index=[time], unit="GtCO2")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC")
    co2_cum = Parameter(index=[time], unit="GtCO2")

    r_0 = Parameter()
    r_C = Parameter()
    r_T = Parameter()
    iIRF_100_Max = Parameter()
    alpha_fit1 = Parameter()
    alpha_fit2 = Parameter()

    function run_timestep(pp, vv, dd, tt)
        vv.iIRF_100[tt] = min(pp.r_0 + pp.r_C * pp.co2_cum[tt] + pp.r_T * pp.T_AT[tt], pp.iIRF_100_Max)
        vv.alpha[tt] = pp.alpha_fit1 * exp(pp.alpha_fit2 * vv.iIRF_100[tt])
    end
end

function addPostTemperature(model, co2calib)
    params = CSV.read("../data/carbon-cycle.csv", DataFrame)

    posttemp = add_comp!(model, PostTemperature)

    posttemp[:r_0] = params[params.Parameter .== "r_{0}", co2calib][1]
    posttemp[:r_C] = params[params.Parameter .== "r_{C}", co2calib][1]
    posttemp[:r_T] = params[params.Parameter .== "r_{T}", co2calib][1]
    posttemp[:iIRF_100_Max] = params[params.Parameter .== "Max iIRF_100", co2calib][1]
    posttemp[:alpha_fit1] = params[params.Parameter .== "alpha fit 1", co2calib][1]
    posttemp[:alpha_fit2] = params[params.Parameter .== "alpha fit 2", co2calib][1]

    posttemp
end
