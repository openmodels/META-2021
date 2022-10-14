@defcomp ISMModel begin
    monsoonsteps = Index()

    # Variables

    BSO4 = Variable(index=[time]) # regional sulphate burden
    Apl = Variable(index=[time]) # planetary albedo
    Aplcrit = Variable(index=[time]) # critical albedo
    mNINO3pt4 = Variable(index=[time], unit="MSLP") # strength of the Walker circulation (MSLP := decadal scale anomaly in mean sea level pressure
    pinit1 = Variable(index=[time]) # rainfall probability in the onset period if Apl < Aplcrit
    pinit = Variable(index=[time]) # rainfall probability in the onset period if Apl >= Aplcrit
    Pwet = Variable(index=[time], unit="mm/day") # daily rainfall intensity
    Pbar = Variable(index=[time], unit="mm/day") # seasonal mean precipitation rate
    D_ISM = Variable(index=[time], unit="pctGDP") # damages as a function of whether the year is drought, flood or neither

    dailyrainfall = Variable(index=[time, monsoonsteps], unit="mm/day")

    # Parameters

    SO_2 = Parameter(index=[time], unit="Mt") # SO2 emissions from the Asia region, from the RCPs
    st_ppm = Parameter(index=[time], unit="ppm") # atmospheric CO2 in ppm, output of CO2 model
    T_AT = Parameter(index=[time], unit="degC") # GMST, output of temperature model
    uniforms = Parameter(index=[time, monsoonsteps])

    Ddrought = Parameter(default=0.035, unit="pctGDP")
    Dflood = Parameter(default=0.0084, unit="pctGDP")
    Pdry = Parameter(default=0, unit="mm/day")
    Pwet_0 = Parameter(default=9, unit="mm/day")
    pinit_0 = Parameter(default=0.75)
    p_m = Parameter(default=0.82)
    ism_delta = Parameter(default=16, unit="days")
    pdoubleprime = Parameter(default=0.42, unit="mm/day")
    p_0 = Parameter(default=0.2)
    pprime = Parameter(default=0.41)
    m_0 = Parameter(default=1008.649, unit="mb")
    mprime = Parameter(default=-0.29)
    mNINO3pt4_0 = Parameter(default=1009.983, unit="mb")
    O = Parameter(default=1.5)
    A_s = Parameter(default=0.16)
    HSO2 = Parameter(default=0.65)
    T_pl = Parameter(default=0.76)
    Beta_pl = Parameter(default=0.29)
    alphapl_3 = Parameter(default=8.5, unit="m^2/g")
    V = Parameter(default=0.005479)
    Omega = Parameter(default=4200000000000, unit="m^2")
    alphapl_1 = Parameter(default=0.02)
    alphapl_2 = Parameter(default=0.37)
    Pdrought = Parameter(default=2.8667, unit="mm/day")
    Pflood = Parameter(default=7.6667, unit="mm/day")

    T_AT_2010 = Parameter(default=0.854, unit="degC")
    f_NINO = Parameter(index=[time])

    function run_timestep(pp, vv, dd, tt)

        vv.BSO4[tt] = pp.SO_2[tt] * 8.5/pp.SO_2[TimestepIndex(1)] * 10 ^ 12 * pp.O * pp.HSO2 * pp.V / pp.Omega
        vv.Apl[tt] = 0.47 + 2 * pp.T_pl ^ 2 * (1 - pp.A_s) ^ 2 * pp.Beta_pl * pp.alphapl_3 * vv.BSO4[tt]
        vv.Aplcrit[tt] = pp.alphapl_1 * log(pp.st_ppm[tt]) + pp.alphapl_2
        vv.mNINO3pt4[tt] = pp.mprime / 30 * (pp.T_AT[tt] - pp.T_AT_2010) * pp.f_NINO[tt] + pp.mNINO3pt4_0
        vv.pinit1[tt] = pp.pprime * (vv.mNINO3pt4[tt] - pp.m_0) + pp.p_0
        vv.pinit[tt] = vv.Apl[tt] < vv.Aplcrit[tt] ? vv.pinit1[tt] : 1 - pp.p_m
        vv.Pwet[tt] = pp.pdoubleprime * (pp.T_AT[tt] - pp.T_AT_2010) + pp.Pwet_0

        # Create a matrix of probabilities of wet/dry days

        if false
            N = 35 # number of periods in monsoon season, each comprising 4 days

            ISM_period1to4_uniform = zeros(4)
            ISM_daily_rainfall_draws = zeros(N)

#            for n in 1:4
#
#                ISM_period1to4_uniform[n] = pp.uniforms[tt,n]
#
#                ISM_period1to4_uniform[n] > 1 - vv.pinit[tt] ? ISM_daily_rainfall_draws[n] = 1 : ISM_daily_rainfall_draws[n] = 0
#
#            end
#
#            for n in 5:N
#
#                ISM_daily_rainfall_draws[n] = pp.uniforms[tt,n]
#
#            end
#
            # Calculate rainfall for each period of days in year t
#            dailyrainfall = zeros(N)

#            for n in 1:4 # rainfall during the onset season of 16 days, i.e., 4 periods

#                ISM_daily_rainfall_draws[n] == 1 ? dailyrainfall[n] = vv.Pwet[tt] : pp.Pdry

#            end

#            for n in 5:N # rainfall during the remainder of the season

#                ISM_daily_rainfall_draws[n] < max(min(1 / pp.ism_delta * (sum(dailyrainfall[n-4:n-1]) * 4 - pp.Pdry) / (vv.Pwet[tt] - pp.Pdry), pp.p_m), 1 - pp.p_m) ? dailyrainfall[n] = vv.Pwet[tt] : pp.Pdry

#            end

#            for n in 1:N
#                vv.dailyrainfall[tt, n] = dailyrainfall[n]
#            end
        else
            for ss in dd.monsoonsteps
                if ss <= 4
                    ISM_daily_rainfall_draws = pp.uniforms[tt, ss] > 1 - vv.pinit[tt]
                else
                    ISM_daily_rainfall_draws = pp.uniforms[tt, ss] < max(min(1 / pp.ism_delta * (sum(vv.dailyrainfall[tt, ss-4:ss-1]) * 4 - pp.Pdry) / (vv.Pwet[tt] - pp.Pdry), pp.p_m), 1 - pp.p_m)
                end

                # Calculate rainfall for each period of days in year t
                vv.dailyrainfall[tt, ss] = ISM_daily_rainfall_draws ? vv.Pwet[tt] : pp.Pdry
            end
        end

        vv.Pbar[tt] = mean(vv.dailyrainfall[tt, :]) # Seasonal mean precipitation rate for year t

        if vv.Pbar[tt] <= pp.Pdrought
            vv.D_ISM[tt] = pp.Ddrought
        elseif vv.Pbar[tt] >= pp.Pflood
            vv.D_ISM[tt] = pp.Dflood
        else
            vv.D_ISM[tt] = 0
        end

    end

end

function addISMModel(model, ismcalib, before=nothing, after=nothing)

    params = CSV.read("../data/ISMparams.csv", DataFrame)

    if ismcalib ∉ names(params)[2:end]
        throw(ArgumentError("Unknown ISM model calibration"))
    end

    ismmodel = add_comp!(model, ISMModel, before=before, after=after)

    ismmodel

end
