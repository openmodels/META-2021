import Mimi.has_parameter

function myupdate_param!(model, comp, param, value)
    unique_name = Symbol("$(comp)_$param")
    if has_parameter(model.md, unique_name)
        update_param!(model, unique_name, value)
    else
        set_param!(model, comp, param, unique_name, value)
    end
end

function preset_fill_notp(model::Model, benchmark::DataFrame, rr::Int64)
    ## Fill in values
    myupdate_param!(model, :CO2Model, :a0, benchmark.a_0[rr])
    myupdate_param!(model, :CO2Model, :a1, benchmark.a_1[rr])
    myupdate_param!(model, :CO2Model, :a3, benchmark.a_3[rr])
    myupdate_param!(model, :CO2Model, :rho1, benchmark.rho_1[rr])
    myupdate_param!(model, :CO2Model, :rho2, benchmark.rho_2[rr])
    myupdate_param!(model, :CO2Model, :rho3, benchmark.rho_3[rr])
    myupdate_param!(model, :PostTemperature, :r_0, benchmark."r_{0) / Distribution"[rr])
    myupdate_param!(model, :PostTemperature, :r_C, benchmark."r_{C} / Distribution"[rr])
    myupdate_param!(model, :PostTemperature, :r_T, benchmark."r_{T} / Distribution"[rr])
    myupdate_param!(model, :TemperatureModel, :xi_1, benchmark.xi_1[rr])
    myupdate_param!(model, :TemperatureModel, :F_2xCO2, benchmark.F_2XCO2[rr])
    myupdate_param!(model, :TemperatureModel, :fair_ECS, benchmark.T_2xCO2[rr])
    myupdate_param!(model, :TemperatureModel, :fair_gamma, benchmark.gamma[rr])
    myupdate_param!(model, :TemperatureModel, :fair_C_0, benchmark.C0[rr])
    myupdate_param!(model, :Forcing, :F_2xCO2, benchmark.F_2XCO2[rr])
    myupdate_param!(model, :CH4Model, :ch4_alpha, benchmark."Dataset 1"[rr])

    myupdate_param!(model, :Consumption, :damagepersist, 0.5)

    slrindexes = findfirst(names(benchmark) .== "Distribution / AFG"):findfirst(names(benchmark) .== "Distribution / ZWE")
    countries = [x[end-2:end] for x in names(benchmark)[slrindexes]]
    slrcoeffvalues = collect(benchmark[rr, slrindexes])
    myupdate_param!(model, :Consumption, :slrcoeff, [slrcoeffvalues[findfirst(countries .== country)] for country in dim_keys(model, :country)])

    myupdate_param!(model, :Utility, :EMUC, 1.5)
    myupdate_param!(model, :Utility, :PRTP, 0.01)
end

function preset_fill_tp(model::Model, benchmark::DataFrame, rr::Int64)
    ## Fill in values
    preset_fill_notp(model, benchmark, rr)

    # Not using distribution
    # myupdate_param!(model, :PCFModel, :propCH4, benchmark."propCH4 / Kessler probabilistic"[rr])
    # myupdate_param!(model, :PCFModel, :beta_PF, benchmark."beta_PF / Kessler probabilistic"[rr])
    # myupdate_param!(model, :PCFModel, :C_PF, benchmark."C_PF (GtC) / Kessler probabilistic"[rr])
    # myupdate_param!(model, :PCFModel, :propPassive, benchmark."propPassive / Kessler probabilistic"[rr])
    # myupdate_param!(model, :PCFModel, :tau, benchmark."tau (years) / Kessler probabilistic"[rr])

    # Not using distribution
    # myupdate_param!(model, :AmazonDieback, :Delta_AMAZ, benchmark."Delta_AMAZ / Distribution"[rr])

    # Not using distribution
    # myupdate_param!(model, :GISModel, :avoldot0, benchmark."avoldot0 / Distribution"[rr])

    myupdate_param!(model, :WAISmodel, :waisrate, benchmark."waisrate / Distribution"[rr] / 1000)

    myupdate_param!(model, :SAFModel, :ECS, benchmark.T_2xCO2[rr])
    if benchmark."Nonlinear SAF / Distribution"[rr] != "Error"
        myupdate_param!(model, :SAFModel, :saf_delta, benchmark."Nonlinear SAF / Distribution"[rr])
    end
    if benchmark."Warming half-life / Distribution"[rr] != "Error"
        myupdate_param!(model, :SAFModel, :FRT, benchmark."Warming half-life / Distribution"[rr])
    end

    include("../src/lib/interacts.jl")

    interactionparams = CSV.read("../data/interactionparams.csv", DataFrame)

    cc1 = findfirst(names(benchmark) .== "Distribution / MGIS")

    function callback(symb, paramaffected, paramforcing)
        row = findfirst((interactionparams."Affected tipping point" .== paramaffected) .& (interactionparams."Forcing tipping point" .== paramforcing))
        myupdate_param!(model, :Interactions, symb, benchmark[rr, cc1 + row - 1])
    end

    allinteractcalls(callback)

    bindrawstarts = findall(x -> occursin("2010 / Binomial draw", x), names(benchmark))
    bindrawends = findall(x -> occursin("2200 / Binomial draw", x), names(benchmark))
    bindrawcomps = [:OMH, :AmazonDieback, :WAISmodel, :AMOC]
    for cc in 1:length(bindrawcomps)
        myupdate_param!(model, bindrawcomps[cc], :uniforms, 1 .- collect(benchmark[rr, bindrawstarts[cc]:bindrawends[cc]]))
    end

    raindrawstart = findfirst(names(benchmark) .== "2010 / Day")
    raindrawend = bindrawstarts[4] - 1
    draws = reshape(collect(benchmark[rr, raindrawstart:raindrawend]), (dim_count(model, :monsoonsteps), dim_count(model, :time)))'
    myupdate_param!(model, :ISMModel, :uniforms, draws)
end
