@defcomp TemperatureModel begin

    # Variables

    # Preliminary temperature estimate
    T_AT1 = Variable(index=[time], unit="degC")
    T_LO1 = Variable(index=[time], unit="degC")
    # Temperature estimate after correction
    T_AT = Variable(index=[time], unit="degC")
    T_LO = Variable(index=[time], unit="degC")

    # Parameters

    F = Parameter(index=[time], unit="W/m^2")
    T_AT_adjustment = Parameter(index=[time], unit="degC")

    xi_1 = Parameter()
    xi_2 = Parameter()
    xi_3 = Parameter()
    xi_4 = Parameter()

    T_AT_2010 = Parameter(unit="degC")
    T_LO_2010 = Parameter(unit="degC")

    function run_timestep(pp, vv, dd, tt)

        if is_first(tt)

            vv.T_AT1[tt] = pp.T_AT_2010
            vv.T_LO1[tt] = pp.T_LO_2010
            vv.T_AT[tt] = pp.T_AT_2010
            vv.T_LO[tt] = pp.T_LO_2010

        else

            vv.T_AT1[tt] = vv.T_AT1[tt-1] + pp.xi_1 * (pp.F[tt] - pp.xi_2 * vv.T_AT1[tt-1] - pp.xi_3 * (vv.T_AT1[tt-1] - vv.T_LO1[tt-1]))
            vv.T_LO1[tt] = vv.T_LO1[tt-1] + pp.xi_4 * (vv.T_AT1[tt-1] - vv.T_LO1[tt-1])

            vv.T_AT[tt] = vv.T_AT1[tt] + pp.T_AT_adjustment[tt]
            vv.T_LO[tt] = vv.T_LO[tt-1] + pp.xi_4 * (vv.T_AT[tt-1] - vv.T_LO[tt-1])

        end

    end

end

function addTemperatureModel(model, tempcalib)

    params = CSV.read("../data/temperatureparams.csv", DataFrame)

    if tempcalib ∉ names(params)[2:end]
        throw(ArgumentError("Unknown temperature model calibration"))
    end

    temperaturemodel = add_comp!(model, TemperatureModel)

    temperaturemodel[:xi_1] = params[params.Parameter .== "xi_1", tempcalib][1]
    temperaturemodel[:xi_2] = params[params.Parameter .== "xi_2", tempcalib][1]
    temperaturemodel[:xi_3] = params[params.Parameter .== "xi_3", tempcalib][1]
    temperaturemodel[:xi_4] = params[params.Parameter .== "xi_4", tempcalib][1]

    temperaturemodel[:T_AT_2010] = params[params.Parameter .== "T_AT_2010", tempcalib][1]
    temperaturemodel[:T_LO_2010] = params[params.Parameter .== "T_LO_2010", tempcalib][1]

    temperaturemodel[:T_AT_adjustment] = zeros(dim_count(model, :time))

    temperaturemodel

end
