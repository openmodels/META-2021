@defcomp Forcing begin

    # Variables

    F = Variable(index=[time], unit="W/m^2")

    # Parameters

    st_ppm = Parameter(index=[time], unit="ppm") # atmospheric CO2 in ppm, output of CO2 model
    F_CH4 = Parameter(index=[time], unit="W/m^2") # forcing from CH4, output of CH4 model [ENSURE CALLED THE SAME]
    F_EX = Parameter(index=[time], unit="W/m^2") # forcing from other agents, from RCPs

    F_2xCO2 = Parameter(unit="W/m^2") # forcing from equilibrium CO2 doubling (Wm^-2)
    ppm_preind = Parameter(unit="ppm") # pre-industrial CO2 concentration (ppm)

    function run_timestep(pp, vv, dd, tt)
        vv.F[tt] = pp.F_2xCO2 * log2(pp.st_ppm[tt]/pp.ppm_preind) + pp.F_CH4[tt] + pp.F_EX[tt]
    end

end

function addForcing(model, tempcalib)

    params = CSV.read("../data/temperatureparams.csv", DataFrame)

    if tempcalib ∉ names(params)[2:end]
        throw(ArgumentError("Unknown temperature model calibration"))
    end

    forcing = add_comp!(model, Forcing)

    forcing[:F_2xCO2] = params[params.Parameter .== "F_2xCO2", tempcalib][1]
    forcing[:ppm_preind] = params[params.Parameter .== "Pre-industrial CO2 (ppm)", tempcalib][1]

    forcing

end
