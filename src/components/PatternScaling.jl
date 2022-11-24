@defcomp PatternScaling begin
    country = Index()

    # Variables
    scale_country = Variable(index=[time, country])
    T_country = Variable(index=[time, country], unit="degC")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC")

    # Country-specific calibration
    ps_alpha = Parameter(index=[country], unit="degC")
    ps_beta = Parameter(index=[country], unit="degC")

    # General parameters
    T_AT_2010 = Parameter(unit="degC", default=0.854) # warming in 2010
    GMST_2010 = Parameter(unit="degC", default=20.02780209) # global mean temperature in 2010

    function run_timestep(pp, vv, dd, tt)
        for cc in dd.country
            if pp.GMST_2010 + pp.T_AT[tt] - pp.T_AT_2010 - 18.825 > 0
                vv.scale_country[tt, cc] = pp.ps_alpha[cc] + pp.ps_beta[cc] * log(pp.GMST_2010 + pp.T_AT[tt] - pp.T_AT_2010 - 18.825)
            else
                vv.scale_country[tt, cc] = pp.ps_alpha[cc] + pp.ps_beta[cc] * log(pp.GMST_2010 - 18.825)
            end
            vv.T_country[tt, cc] = (pp.GMST_2010+pp.T_AT[tt]-pp.T_AT[TimestepIndex(1)])*vv.scale_country[tt, cc]
        end
    end
end

function addPatternScaling(model)
    params = CSV.read("../data/pattern-scaling.csv", DataFrame)

    pattscale = add_comp!(model, PatternScaling)

    pattscale[:ps_alpha] = params.alpha
    pattscale[:ps_beta] = params.beta

    pattscale
end
