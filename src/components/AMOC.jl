@defcomp AMOC begin
    country = Index()

    # Variables
    p_AMOC = Variable(index=[time], unit="fraction")
    I_AMOC = Variable{Bool}(index=[time])
    deltaT_country_AMOC = Variable(index=[time, country], unit="degC")
    T_country_AMOC = Variable(index=[time, country], unit="degC")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC")
    f_AMOC = Parameter(index=[time])

    uniforms = Parameter(index=[time])
    b_AMOC = Parameter()
    Delta_AMOC = Parameter(unit="year", default=35)
    GMST_2010 = Parameter(unit="degC", default=20.02780209) # Source: RCP4.5
    max_deltaT_country_AMOC = Parameter(index=[country], unit="degC")

    ## Additional parameters for computing T_country_AMOC
    scale_country = Parameter(index=[time, country], unit="degC")

    function run_timestep(pp, vv, dd, tt)
        vv.p_AMOC[tt] = min((1-exp(-pp.b_AMOC*pp.T_AT[tt]))*pp.f_AMOC[tt], 1)

        if is_first(tt)
            vv.I_AMOC[tt] = pp.uniforms[tt] < vv.p_AMOC[tt]

            for cc in dd.country
                vv.deltaT_country_AMOC[tt, cc] = vv.I_AMOC[tt] ? pp.max_deltaT_country_AMOC[cc] / pp.Delta_AMOC : 0
            end
        else
            vv.I_AMOC[tt] = vv.I_AMOC[tt-1] || (pp.uniforms[tt] < vv.p_AMOC[tt])

            for cc in dd.country
                if vv.I_AMOC[tt] && abs(vv.deltaT_country_AMOC[tt-1, cc]) < abs(pp.max_deltaT_country_AMOC[cc])
                    vv.deltaT_country_AMOC[tt, cc] = vv.deltaT_country_AMOC[tt-1, cc] + pp.max_deltaT_country_AMOC[cc] / pp.Delta_AMOC
                    if abs(vv.deltaT_country_AMOC[tt, cc]) > abs(pp.max_deltaT_country_AMOC[cc])
                        vv.deltaT_country_AMOC[tt, cc] = pp.max_deltaT_country_AMOC[cc]
                    end
                else
                    vv.deltaT_country_AMOC[tt, cc] = vv.deltaT_country_AMOC[tt-1, cc]
                end
            end
        end

        for cc in dd.country
            vv.T_country_AMOC[tt, cc] = (pp.GMST_2010+pp.T_AT[tt]-pp.T_AT[TimestepIndex(1)])*pp.scale_country[tt, cc] + vv.deltaT_country_AMOC[tt, cc]
        end
    end
end

b_AMOC_calibs = Dict{String, Float64}("Hadley" => 0.135, "BCM" => 0.611,
	                              "IPSL" => 0.54, "HADCM" => 1.6)

function addAMOC(model, calibration; before=nothing, after=nothing)
    if calibration ∉ keys(b_AMOC_calibs)
        throw(ArgumentError("Unknown AMOC model calibration"))
    end

    params = CSV.read("../data/AMOCparams.csv", DataFrame)

    amoc = add_comp!(model, AMOC, before=before, after=after)
    amoc[:b_AMOC] = b_AMOC_calibs[calibration]
    amoc[:max_deltaT_country_AMOC] = params[!, calibration]

    amoc
end
