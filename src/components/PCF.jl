@defcomp PCFModel begin

    # Variables

    PF_extent = Variable(index=[time], unit="percent") # area of permafrost remaining relative to time zero
    C_thawedPF = Variable(index=[time], unit="GtC") # amount of carbon in freshly thawed permafrost (GtC)
    CCum_PF = Variable(index=[time], unit="GtC") # cumulative CO2 emissions from thawed permafrost (GtC)
    DifferencedCO2 = Variable(index=[time], unit="GtC") # difference in cumulative permafrost CO2 emissions (GtC)
    CO2_PF = Variable(index=[time], unit="GtCO2") # permafrost CO2 emissions (GtCO2)
    CH4_PF = Variable(index=[time], unit="MtCH4") # permafrost CH4 emissions (MtCH4)

    # Parameters

    T_AT = Parameter(index=[time], unit="degC") # Atmospheric temperature relative to pre-industrial, output of temperature model

    beta_PF = Parameter() # coefficient representing the sensitivity of permafrost thaw to temperature
    C_PF = Parameter(unit="GtC") # total stock of carbon locked in the near-surface northern circumpolar permafrost region
    propPassive = Parameter(unit="percent") # proportion of thawed permafrost in the passive reservoir
    tau = Parameter(unit="years") # e-folding time of permafrost decomposition in the active reservoir
    propCH4 = Parameter(unit="percent") # share of CH4 emissions in total carbon emissions

    T_AT_2010 = Parameter(default=0.854, unit="degC")

    function run_timestep(pp, vv, dd, tt)

        vv.PF_extent[tt] = 1 - pp.beta_PF * (pp.T_AT[tt] - pp.T_AT_2010)

        if is_first(tt)

            vv.C_thawedPF[tt] = 0
            vv.CCum_PF[tt] = 0
            vv.CO2_PF[tt] = 0
            vv.CH4_PF[tt] = 0

        else

            vv.C_thawedPF[tt] = -pp.C_PF * (vv.PF_extent[tt] - vv.PF_extent[tt-1])
            DecompPF = [vv.C_thawedPF[tt-dt] * (1 - pp.propPassive) * (1 - exp(-dt / pp.tau)) for dt in 0:tt.t-1]
            vv.CCum_PF[tt] = sum(DecompPF)
            vv.DifferencedCO2[tt] = vv.CCum_PF[tt] - vv.CCum_PF[tt-1]
            vv.CO2_PF[tt] = vv.DifferencedCO2[tt] * (1 - pp.propCH4) * 44/12
            vv.CH4_PF[tt] = vv.DifferencedCO2[tt] * pp.propCH4 * 16 / 12 * 1000

        end


    end
end

function addPCFModel(model, pcfcalib; before=nothing, after=nothing)

    params = CSV.read("../data/PCFparams.csv", DataFrame)

    if pcfcalib ∉ names(params)[2:end-1]
        throw(ArgumentError("Unknown pcf model calibration"))
    end

    pcfmodel = add_comp!(model, PCFModel, before=before, after=after)

    pcfmodel[:beta_PF] = params[params.Parameter .== "beta_PF", pcfcalib][1]
    pcfmodel[:C_PF] = params[params.Parameter .== "C_PF", pcfcalib][1]
    pcfmodel[:propPassive] = params[params.Parameter .== "propPassive", pcfcalib][1]
    pcfmodel[:tau] = params[params.Parameter .== "tau", pcfcalib][1]
    pcfmodel[:propCH4] = params[params.Parameter .== "propCH4", pcfcalib][1]

    # pcfmodel[:T_AT_2010] = params[params.Parameter .== "T_AT_2010", pcfcalib][1]

    pcfmodel

end
