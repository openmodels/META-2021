@defcomp CO2Model begin
    # Variables
    s0 = Variable(index=[time], unit="GtCO2")
    s1 = Variable(index=[time], unit="GtCO2")
    s2 = Variable(index=[time], unit="GtCO2")
    s3 = Variable(index=[time], unit="GtCO2")
    st = Variable(index=[time], unit="GtCO2")
    st_ppm = Variable(index=[time], unit="ppm")
    co2_cum_ppm = Variable(index=[time], unit="ppm")
    co2_cum = Variable(index=[time], unit="GtCO2")
    co2_total = Variable(index=[time], unit="GtCO2")

    # Parameters
    co2_2009 = Parameter(unit="GtCO2", default=35.70037753)
    co2_rcp = Parameter(index=[time], unit="GtCO2")
    co2_pcf = Parameter(index=[time], unit="GtCO2")
    co2_amazon = Parameter(index=[time], unit="GtCO2")
    co2_extra = Parameter(index=[time], unit="GtCO2")

    alpha = Parameter(index=[time], unit="GtCO2")

    a0 = Parameter()
    a1 = Parameter()
    a2 = Parameter()
    a3 = Parameter()
    rho0 = Parameter()
    rho1 = Parameter()
    rho2 = Parameter()
    rho3 = Parameter()

    s0_ppm_2010 = Parameter(unit="ppm")
    s1_ppm_2010 = Parameter(unit="ppm")
    s2_ppm_2010 = Parameter(unit="ppm")
    s3_ppm_2010 = Parameter(unit="ppm")
    ppm_preind = Parameter(unit="ppm")
    ppm_to_gtco2 = Parameter(unit="GtCO2/ppm")
    c_to2010 = Parameter(unit="GtC")

    function run_timestep(pp, vv, dd, tt)
        if is_first(tt)
            vv.s0[tt] = pp.s0_ppm_2010 * pp.ppm_to_gtco2
            vv.s1[tt] = pp.s1_ppm_2010 * pp.ppm_to_gtco2
            vv.s2[tt] = pp.s2_ppm_2010 * pp.ppm_to_gtco2
            vv.s3[tt] = pp.s3_ppm_2010 * pp.ppm_to_gtco2

            vv.co2_total[tt] = pp.co2_2009
            vv.co2_cum_ppm[tt] = vv.co2_total[tt]
        else
            vv.co2_total[tt] = pp.co2_rcp[tt-1] / 1000 + pp.co2_pcf[tt-1] + pp.co2_amazon[tt-1] + pp.co2_extra[tt]

            vv.s0[tt] = pp.a0 * vv.co2_total[tt] + (1 - pp.rho0/pp.alpha[tt-1]) * vv.s0[tt-1]
            vv.s1[tt] = pp.a1 * vv.co2_total[tt] + (1 - pp.rho1/pp.alpha[tt-1]) * vv.s1[tt-1]
            vv.s2[tt] = pp.a2 * vv.co2_total[tt] + (1 - pp.rho2/pp.alpha[tt-1]) * vv.s2[tt-1]
            vv.s3[tt] = pp.a3 * vv.co2_total[tt] + (1 - pp.rho3/pp.alpha[tt-1]) * vv.s3[tt-1]

            vv.co2_cum_ppm[tt] = vv.co2_total[tt] + vv.co2_cum_ppm[tt-1]
        end

        vv.st[tt] = vv.s0[tt] + vv.s1[tt] + vv.s2[tt] + vv.s3[tt] + pp.ppm_preind * pp.ppm_to_gtco2
        vv.st_ppm[tt] = vv.st[tt] / pp.ppm_to_gtco2
        vv.co2_cum[tt] = vv.co2_cum_ppm[tt] * 12/44 + pp.c_to2010 - (vv.s0[tt] + vv.s1[tt] + vv.s2[tt] + vv.s3[tt]) * 12/44
    end
end

function addCO2Model(model, co2calib)
    if co2calib == "Distribution"
        error("Distribution CO2 cycle not implemented")
    end

    params = CSV.read("../data/carbon-cycle.csv", DataFrame)
    if co2calib ∉ names(params)[2:end-1]
        throw(ArgumentError("Unknown CO2 cycle calibration"))
    end

    co2model = add_comp!(model, CO2Model)

    co2model[:co2_pcf] = zeros(dim_count(model, :time))
    co2model[:co2_amazon] = zeros(dim_count(model, :time))
    co2model[:co2_extra] = zeros(dim_count(model, :time))

    co2model[:a0] = params[params.Parameter .== "a_0", co2calib][1]
    co2model[:a1] = params[params.Parameter .== "a_1", co2calib][1]
    co2model[:a2] = params[params.Parameter .== "a_2", co2calib][1]
    co2model[:a3] = params[params.Parameter .== "a_3", co2calib][1]
    co2model[:rho0] = params[params.Parameter .== "rho_0", co2calib][1]
    co2model[:rho1] = params[params.Parameter .== "rho_1", co2calib][1]
    co2model[:rho2] = params[params.Parameter .== "rho_2", co2calib][1]
    co2model[:rho3] = params[params.Parameter .== "rho_3", co2calib][1]

    co2model[:s0_ppm_2010] = params[params.Parameter .== "S_0 in 2010 (ppm)", co2calib][1]
    co2model[:s1_ppm_2010] = params[params.Parameter .== "S_1 in 2010 (ppm)", co2calib][1]
    co2model[:s2_ppm_2010] = params[params.Parameter .== "S_2 in 2010 (ppm)", co2calib][1]
    co2model[:s3_ppm_2010] = params[params.Parameter .== "S_3 in 2010 (ppm)", co2calib][1]
    co2model[:ppm_preind] = params[params.Parameter .== "Pre-industrial CO2 (ppm)", co2calib][1]
    co2model[:ppm_to_gtco2] = params[params.Parameter .== "Conversion factor ppm to GtCO2", co2calib][1]
    co2model[:c_to2010] = params[params.Parameter .== "Cumulative CO2 emissions, 1750-2010 (GtC)", co2calib][1]

    co2model
end
