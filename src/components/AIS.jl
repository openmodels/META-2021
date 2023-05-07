@defcomp AISmodel begin

    # Variables

    # Surface Mass Balance module

    ΔA = Variable(index=[time], unit="percent") # change in continental-scale accumulation from 2010 (%)
    ΔSMB = Variable(index=[time], unit="m") # unadjusted annual mass change (m)
    Adjustment = Variable(index=[time], unit="m") # adjustment factor (m)
    AdjustedΔSMB = Variable(index=[time], unit="m") # adjusted annual mass change (m)
    totalAdjustedΔSMB = Variable(index=[time], unit="m") # cumulative adjusted annual mass change (m)

    # Ice sheet module

    ΔT₀_EAIS = Variable(index=[time], unit="degC") # subsurface oceanic warming in EAIS region of Antarctica (K)
    ΔM_EAIS = Variable(index=[time], unit="m/year") # basal ice shelf melting rate in EAIS region of Antarctica (m/year)
    ΔT₀_Ross = Variable(index=[time], unit="deg")
    ΔM_Ross = Variable(index=[time], unit="m/year")
    ΔT₀_Amundsen = Variable(index=[time], unit="degC")
    ΔM_Amundsen = Variable(index=[time], unit="m/year")
    ΔT₀_Weddell = Variable(index=[time], unit="degC")
    ΔM_Weddell = Variable(index=[time], unit="m/year")
    ΔT₀_Peninsula = Variable(index=[time], unit="degC")
    ΔM_Peninsula = Variable(index=[time], unit="m/year")
    marginalSLR_EAIS = Variable(index=[time], unit="m/year") # marginal SLR in EAIS region of Antarctica (m/year)
    totalSLR_EAIS = Variable(index=[time], unit="m") # total SLR in EAIS region of Antarctica (m)
    marginalSLR_Ross = Variable(index=[time], unit="m/year")
    totalSLR_Ross = Variable(index=[time], unit="m")
    marginalSLR_Amundsen = Variable(index=[time], unit="m/year")
    totalSLR_Amundsen = Variable(index=[time], unit="m")
    marginalSLR_Weddell = Variable(index=[time], unit="m/year")
    totalSLR_Weddell = Variable(index=[time], unit="m")
    marginalSLR_Peninsula = Variable(index=[time], unit="m/year")
    totalSLR_Peninsula = Variable(index=[time], unit="m")

    SLR_AIS = Variable(index=[time], unit="m") # total SLR from all AIS sources (m)

    # Parameters

    T_AT = Parameter(index=[time], unit="degC") # GMST increase (degC), series starting in 2010
    T_AT_tminus100 = Parameter(index=[time], unit="degC") # GMST increase (degC), series starting in 1910
    β_EAIS = Parameter(unit="degC/degC") # scaling coefficient
    δ_EAIS = Parameter{Int64}(unit="years") # time delay (years)
    β_Ross = Parameter(unit="degC/degC")
    δ_Ross = Parameter{Int64}(unit="years")
    β_Amundsen = Parameter(unit="degC/degC")
    δ_Amundsen = Parameter{Int64}(unit="years")
    β_Weddell = Parameter(unit="degC/degC")
    δ_Weddell = Parameter{Int64}(unit="years")
    β_Peninsula = Parameter(unit="degC/degC")
    δ_Peninsula = Parameter{Int64}(unit="years")
    R_functions_EAIS = Parameter(index=[time]) # SLR response function values
    R_functions_Ross = Parameter(index=[time])
    R_functions_Amundsen = Parameter(index=[time])
    R_functions_Weddell = Parameter(index=[time])
    R_functions_Peninsula = Parameter(index=[time])

    # Surface Mass Balance module

    ϕ = Parameter(default=1.2) # scaling coefficient global to regional temperature
    ω = Parameter(default=-0.05) # change in continental accumulation per degree of Antarctic warming (% K^{-1})
    γ = Parameter(default=0.00795) # interaction between SMB and dynamic processes
    K = Parameter(default=0.008) # SMB eventually after threshold crossed (m yr^{-1})
    C = Parameter(default=1) # parameter of generalized logistic function
    Q = Parameter(default=0.5) # parameter of generalized logistic function
    B = Parameter(default=1.5) # parameter of generalized logistic function
    V = Parameter(default=0.5) # parameter of generalized logistic function

    # Ice sheet module

    λ = Parameter(default=11.5) # basal melt sensitivity parameter (ma^{-1} K^{-1})

    function run_timestep(pp, vv, dd, tt)
        if is_first(tt)

            vv.totalAdjustedΔSMB[tt] = 0

            vv.ΔT₀_EAIS[tt] = (gettime(tt) - pp.δ_EAIS < 2010 ? pp.β_EAIS * pp.T_AT_tminus100[tt + 100 - pp.δ_EAIS] : pp.β_EAIS * pp.T_AT[tt - pp.δ_EAIS])
            vv.ΔM_EAIS[tt] = pp.λ * vv.ΔT₀_EAIS[tt]
            vv.marginalSLR_EAIS[tt] = pp.R_functions_EAIS[tt] * vv.ΔM_EAIS[tt]
            vv.totalSLR_EAIS[tt] = vv.marginalSLR_EAIS[tt]

            vv.ΔT₀_Ross[tt] = (gettime(tt) - pp.δ_Ross < 2010 ? pp.β_Ross * pp.T_AT_tminus100[tt + 100 - pp.δ_Ross] : pp.β_Ross * pp.T_AT[tt - pp.δ_Ross])
            vv.ΔM_Ross[tt] = pp.λ * vv.ΔT₀_Ross[tt]
            vv.marginalSLR_Ross[tt] = pp.R_functions_Ross[tt] * vv.ΔM_Ross[tt]
            vv.totalSLR_Ross[tt] = vv.marginalSLR_Ross[tt]

            vv.ΔT₀_Amundsen[tt] = (gettime(tt) - pp.δ_Amundsen < 2010 ? pp.β_Amundsen * pp.T_AT_tminus100[tt + 100 - pp.δ_Amundsen] : pp.β_Amundsen * pp.T_AT[tt - pp.δ_Amundsen])
            vv.ΔM_Amundsen[tt] = pp.λ * vv.ΔT₀_Amundsen[tt]
            vv.marginalSLR_Amundsen[tt] = pp.R_functions_Amundsen[tt] * vv.ΔM_Amundsen[tt]
            vv.totalSLR_Amundsen[tt] = vv.marginalSLR_Amundsen[tt]

            vv.ΔT₀_Weddell[tt] = (gettime(tt) - pp.δ_Weddell < 2010 ? pp.β_Weddell * pp.T_AT_tminus100[tt + 100 - pp.δ_Weddell] : pp.β_Weddell * pp.T_AT[tt - pp.δ_Weddell])
            vv.ΔM_Weddell[tt] = pp.λ * vv.ΔT₀_Weddell[tt]
            vv.marginalSLR_Weddell[tt] = pp.R_functions_Weddell[tt] * vv.ΔM_Weddell[tt]
            vv.totalSLR_Weddell[tt] = vv.marginalSLR_Weddell[tt]

            vv.ΔT₀_Peninsula[tt] = (gettime(tt) - pp.δ_Peninsula < 2010 ? pp.β_Peninsula * pp.T_AT_tminus100[tt + 100 - pp.δ_Peninsula] : pp.β_Peninsula * pp.T_AT[tt - pp.δ_Peninsula])
            vv.ΔM_Peninsula[tt] = pp.λ * vv.ΔT₀_Peninsula[tt]
            vv.marginalSLR_Peninsula[tt] = pp.R_functions_Peninsula[tt] * vv.ΔM_Peninsula[tt]
            vv.totalSLR_Peninsula[tt] = vv.marginalSLR_Peninsula[tt]

        else

            vv.ΔA[tt] = pp.ϕ * pp.ω * ( pp.T_AT[tt] - pp.T_AT[TimestepIndex(1)] )
            vv.ΔSMB[tt] = pp.γ * ( gettime(tt) - 2010 ) ^ -0.1 * vv.ΔA[tt]
            vv.Adjustment[tt] = ( pp.K - vv.ΔSMB[tt] ) / ( pp.C + pp.Q * exp( -pp.B * ( pp.T_AT[tt] - 6.75 ) ) ) ^ ( 1 / pp.V )
            vv.AdjustedΔSMB[tt] = vv.ΔSMB[tt] + vv.Adjustment[tt]
            vv.totalAdjustedΔSMB[tt] = vv.totalAdjustedΔSMB[tt-1] + vv.AdjustedΔSMB[tt]

            vv.ΔT₀_EAIS[tt] = (gettime(tt) - pp.δ_EAIS < 2010 ? pp.β_EAIS * pp.T_AT_tminus100[tt + 100 - pp.δ_EAIS] : pp.β_EAIS * pp.T_AT[tt - pp.δ_EAIS])
            vv.ΔM_EAIS[tt] = pp.λ * vv.ΔT₀_EAIS[tt]
            vv.marginalSLR_EAIS[tt] = pp.R_functions_EAIS[tt] * vv.ΔM_EAIS[tt]
            vv.totalSLR_EAIS[tt] = vv.totalSLR_EAIS[tt-1] + vv.marginalSLR_EAIS[tt]

            vv.ΔT₀_Ross[tt] = (gettime(tt) - pp.δ_Ross < 2010 ? pp.β_Ross * pp.T_AT_tminus100[tt + 100 - pp.δ_Ross] : pp.β_Ross * pp.T_AT[tt - pp.δ_Ross])
            vv.ΔM_Ross[tt] = pp.λ * vv.ΔT₀_Ross[tt]
            vv.marginalSLR_Ross[tt] = pp.R_functions_Ross[tt] * vv.ΔM_Ross[tt]
            vv.totalSLR_Ross[tt] = vv.totalSLR_Ross[tt-1] + vv.marginalSLR_Ross[tt]

            vv.ΔT₀_Amundsen[tt] = (gettime(tt) - pp.δ_Amundsen < 2010 ? pp.β_Amundsen * pp.T_AT_tminus100[tt + 100 - pp.δ_Amundsen] : pp.β_Amundsen * pp.T_AT[tt - pp.δ_Amundsen])
            vv.ΔM_Amundsen[tt] = pp.λ * vv.ΔT₀_Amundsen[tt]
            vv.marginalSLR_Amundsen[tt] = pp.R_functions_Amundsen[tt] * vv.ΔM_Amundsen[tt]
            vv.totalSLR_Amundsen[tt] = vv.totalSLR_Amundsen[tt-1] + vv.marginalSLR_Amundsen[tt]

            vv.ΔT₀_Weddell[tt] = (gettime(tt) - pp.δ_Weddell < 2010 ? pp.β_Weddell * pp.T_AT_tminus100[tt + 100 - pp.δ_Weddell] : pp.β_Weddell * pp.T_AT[tt - pp.δ_Weddell])
            vv.ΔM_Weddell[tt] = pp.λ * vv.ΔT₀_Weddell[tt]
            vv.marginalSLR_Weddell[tt] = pp.R_functions_Weddell[tt] * vv.ΔM_Weddell[tt]
            vv.totalSLR_Weddell[tt] = vv.totalSLR_Weddell[tt-1] + vv.marginalSLR_Weddell[tt]

            vv.ΔT₀_Peninsula[tt] = (gettime(tt) - pp.δ_Peninsula < 2010 ? pp.β_Peninsula * pp.T_AT_tminus100[tt + 100 - pp.δ_Peninsula] : pp.β_Peninsula * pp.T_AT[tt - pp.δ_Peninsula])
            vv.ΔM_Peninsula[tt] = pp.λ * vv.ΔT₀_Peninsula[tt]
            vv.marginalSLR_Peninsula[tt] = pp.R_functions_Peninsula[tt] * vv.ΔM_Peninsula[tt]
            vv.totalSLR_Peninsula[tt] = vv.totalSLR_Peninsula[tt-1] + vv.marginalSLR_Peninsula[tt]

        end

        vv.SLR_AIS[tt] = vv.totalAdjustedΔSMB[tt] + vv.totalSLR_EAIS[tt] + vv.totalSLR_Ross[tt] + vv.totalSLR_Amundsen[tt] + vv.totalSLR_Weddell[tt] + vv.totalSLR_Peninsula[tt]

    end

end

function addAISmodel(model; before=nothing, after=nothing)
    aisgcms = CSV.read("../data/Basal_melt_models.csv", DataFrame)
    aisresponse_EAIS = CSV.read("../data/Response functions - EAIS.csv", DataFrame)
    aisresponse_Ross = CSV.read("../data/Response functions - Ross.csv", DataFrame)
    aisresponse_Amundsen = CSV.read("../data/Response functions - Amundsen.csv", DataFrame)
    aisresponse_Weddell = CSV.read("../data/Response functions - Weddell.csv", DataFrame)
    aisresponse_Peninsula = CSV.read("../data/Response functions - Peninsula.csv", DataFrame)

    ais = add_comp!(model, AISmodel, first=2010, before=before, after=after)

    ais[:β_EAIS] = aisgcms.EAIS_beta[16]
    ais[:δ_EAIS] = aisgcms.EAIS_delta[16]
    ais[:β_Ross] = aisgcms.Ross_beta[14]
    ais[:δ_Ross] = aisgcms.Ross_delta[7]
    ais[:β_Amundsen] = aisgcms.Amundsen_beta[15]
    ais[:δ_Amundsen] = aisgcms.Amundsen_delta[5]
    ais[:β_Weddell] = aisgcms.Weddell_beta[17]
    ais[:δ_Weddell] = aisgcms.Weddell_delta[17]
    ais[:β_Peninsula] = aisgcms.Peninsula_beta[17]
    ais[:δ_Peninsula] = aisgcms.Peninsula_delta[17]
    ais[:R_functions_EAIS] = [zeros(sum(dim_keys(model, :time) .< 2010)); aisresponse_EAIS[!, 2 + 1]]
    ais[:R_functions_Ross] = [zeros(sum(dim_keys(model, :time) .< 2010)); aisresponse_Ross[!, 14 + 1]]
    ais[:R_functions_Amundsen] = [zeros(sum(dim_keys(model, :time) .< 2010)); aisresponse_Amundsen[!, 7 + 1]]
    ais[:R_functions_Weddell] = [zeros(sum(dim_keys(model, :time) .< 2010)); aisresponse_Weddell[!, 8 + 1]]
    ais[:R_functions_Peninsula] = [zeros(sum(dim_keys(model, :time) .< 2010)); aisresponse_Peninsula[!, 11 + 1]]

    ais
end
