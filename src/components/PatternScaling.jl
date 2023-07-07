@defcomp PatternScaling begin
    country = Index()

    # Variables
    # scale_country = Variable(index=[time, country])
    T_country = Variable(index=[time, country], unit="degC")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC")

    # Country-specific calibration
    ps_alpha = Parameter(index=[country], unit="degC")
    ps_beta = Parameter(index=[country], unit="degC")

    function run_timestep(pp, vv, dd, tt)
        for cc in dd.country
            vv.T_country[tt, cc] = pp.ps_alpha[cc] + pp.ps_beta[cc] * pp.T_AT[tt]
        end
    end
end

function addPatternScaling(model)
    params = CSV.read("../data/pattern-scaling_new.csv", DataFrame)

    pattscale = add_comp!(model, PatternScaling, first=2010)

    pattscale[:ps_alpha] = params.alpha
    pattscale[:ps_beta] = params.beta

    pattscale
end
