@defcomp SAFModel begin

    # Variables

    T_AT_PAGE1 = Variable(index=[time], unit="degC") # PAGE-ICE GMST
    SAF_approx = Variable(index=[time], unit="W/m^2") # surface albedo feedback using a quadratic approximation
    ECS_adjusted = Variable(index=[time], unit="degC") # adjusted ECS
    FRT_adjusted = Variable(index=[time], unit="years") # adjusted warming halflife
    T_AT_PAGE2 = Variable(index=[time], unit="degC") # adjusted temperature time-series
    deltaFSAF = Variable(index=[time], unit="W/m^2") # adjustment to SAF forcing
    T_AT_adjustment = Variable(index=[time], unit="degC") # adjustment to GMST

    # Parameters

    F = Parameter(index=[time], unit="W/m^2")

    SAF_bar = Parameter(default=0.349590471, unit="W/m^2/degC") # constant approximation to the SAF
    beta_2 = Parameter(default=-0.001958477, unit="W/m^2/K^3") # T^2 coefficient for the SAF quadratic
    beta_1 = Parameter(default=-0.006387192, unit="W/m^2/K^2") # T^1 coefficient for the SAF quadratic
    beta_0 = Parameter(default=0.363696355, unit="W/m^2/K") # T^0 coefficient for the SAF quadratic
    gamma = Parameter(default=0.109096643, unit="W/m^2/K") # standard deviation of the SAF quadratic
    F_sl = Parameter(default=5.5, unit="W/m^2") # forcing slope
    saf_delta = Parameter() # non-linearity of the SAF
    FRT = Parameter(unit="years") # warming halflife
    saf_alpha = Parameter(default=0.063605129) # linear SAF segment mean
    sigma = Parameter(default=0.029077492) # linear SAF segment standard deviation
    ECS = Parameter()

    T_AT_2010 = Parameter(default=0.854, unit="degC")

    function run_timestep(pp, vv, dd, tt)

        psi = pp.beta_2 * (10 ^ 3) / 3 + pp.beta_1 * (10 ^ 2) / 2 + pp.beta_0 * 10 + pp.gamma * 10 * pp.saf_delta # integration constant for SAF forcing at the segment switch point
        FSAF_0 = pp.beta_2 * (pp.T_AT_2010 ^ 3) / 3 + pp.beta_1 * (pp.T_AT_2010 ^ 2) / 2 + pp.beta_0 * pp.T_AT_2010 + pp.gamma * pp.T_AT_2010 * pp.saf_delta # base year SAF forcing (W/m2)

        if is_first(tt)

            vv.T_AT_PAGE1[tt] = pp.T_AT_2010
            vv.ECS_adjusted[tt] = pp.ECS
            vv.T_AT_PAGE2[tt] = pp.T_AT_2010

        else
            vv.T_AT_PAGE1[tt] = vv.T_AT_PAGE1[tt-1] + (pp.ECS * pp.F[tt-1] / (pp.F_sl * log(2)) - pp.FRT * pp.ECS * (pp.F[tt] - pp.F[tt-1]) / (pp.F_sl * log(2)) - vv.T_AT_PAGE1[tt-1]) * (1 - exp(-1 / pp.FRT)) + pp.ECS * (pp.F[tt] - pp.F[tt-1]) / (pp.F_sl * log(2))
            vv.SAF_approx[tt] = ((pp.beta_2 * (vv.T_AT_PAGE1[tt] ^ 3) / 3 + pp.beta_1 * (vv.T_AT_PAGE1[tt] ^ 2) / 2 + pp.beta_0 * vv.T_AT_PAGE1[tt]) + pp.gamma * vv.T_AT_PAGE1[tt] * pp.saf_delta - FSAF_0) / (vv.T_AT_PAGE1[tt] - pp.T_AT_2010)
            vv.ECS_adjusted[tt] = pp.ECS / (1 - pp.ECS * (vv.SAF_approx[tt] - pp.SAF_bar) / (log(2) * pp.F_sl))
            vv.FRT_adjusted[tt] = pp.FRT / (1 - pp.ECS * (vv.SAF_approx[tt] - pp.SAF_bar) / (log(2) * pp.F_sl))

            if vv.T_AT_PAGE2[tt-1] < 10

                vv.deltaFSAF[tt] = pp.beta_2 * (vv.T_AT_PAGE2[tt-1] ^ 3) / 3 + pp.beta_1 * (vv.T_AT_PAGE2[tt-1] ^ 2) / 2 + pp.beta_0 * vv.T_AT_PAGE2[tt-1] + pp.gamma * vv.T_AT_PAGE2[tt-1] * pp.saf_delta - vv.SAF_approx[tt] * vv.T_AT_PAGE2[tt-1]

            else

                vv.deltaFSAF[tt] = psi + pp.saf_alpha * (vv.T_AT_PAGE2[tt-1] - 10) + pp.sigma * (vv.T_AT_PAGE2[tt-1] - 10) * pp.saf_delta - vv.SAF_approx[tt] * vv.T_AT_PAGE2[tt-1]

            end

            vv.T_AT_PAGE2[tt] = vv.T_AT_PAGE2[tt-1] + (vv.ECS_adjusted[tt] * (pp.F[tt-1] - (pp.F[tt] - pp.F[tt-1]) * vv.FRT_adjusted[tt] + vv.deltaFSAF[tt]) / (log(2) * pp.F_sl) - vv.T_AT_PAGE2[tt-1]) * (1 - exp(-1 / vv.FRT_adjusted[tt])) + vv.ECS_adjusted[tt] * (pp.F[tt] - pp.F[tt-1]) / (log(2) * pp.F_sl)

        end

       vv.T_AT_adjustment[tt] = vv.T_AT_PAGE2[tt] - vv.T_AT_PAGE1[tt]

    end

end

function addSAFModel(model, safcalib; before=nothing, after=nothing)

    params = CSV.read("../data/SAFparams.csv", DataFrame)

    if safcalib ∉ names(params)[2:end]
        throw(ArgumentError("Unknown SAF model calibration"))
    end

    safmodel = add_comp!(model, SAFModel, before=before, after=after)

    #safmodel[:SAF_bar] = params[params.Parameter .== "ECS-average SAF (W/m2/degC)", safcalib][1]
    #safmodel[:beta_2] = params[params.Parameter .== "SAF_mean_quad_segment: T^2 coeff (W/m2/K3)", safcalib][1]
    #safmodel[:beta_1] = params[params.Parameter .== "SAF_mean_quad_segment: T^1 coeff (W/m2/K2)", safcalib][1]
    #safmodel[:beta_0] = params[params.Parameter .== "SAF_mean_quad_segment: T^0 coeff (W/m2/K)", safcalib][1]
    #safmodel[:gamma] = params[params.Parameter .== "SAF_std_quad_segment (W/m2/K)", safcalib][1]
    #safmodel[:F_sl] = params[params.Parameter .== "Forcing slope (W/m2)", safcalib][1]
    safmodel[:saf_delta] = params[params.Parameter .== "Nonlinear SAF", safcalib][1]
    safmodel[:FRT] = params[params.Parameter .== "Warming half-life", safcalib][1]
    #safmodel[:saf_alpha] = params[params.Parameter .== "Linear SAF Mean", safcalib][1]
    #safmodel[:sigma] = params[params.Parameter .== "Linear SAF Std. Dev.", safcalib][1]
    safmodel[:ECS] = params[params.Parameter .== "ECS", safcalib][1]

    #safmodel[:T_AT_2010] = params[params.Parameter .== "T_AT_2010", safcalib][1]

    safmodel

end
