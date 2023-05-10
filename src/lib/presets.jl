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
    beta1indexes = findall(x -> occursin("beta1dist", x), names(benchmark))
    countries = [x[1:3] for x in names(benchmark)[beta1indexes]]
    beta1s = collect(benchmark[rr, beta1indexes])
    beta2indexes = findall(x -> occursin("beta2dist", x), names(benchmark))
    beta2s = collect(benchmark[rr, beta2indexes])

    myupdate_param!(model, :Consumption, :seeds, zeros(dim_count(model, :country)))
    myupdate_param!(model, :Consumption, :beta1, [beta1s[findfirst(countries .== country)] for country in dim_keys(model, :country)])
    myupdate_param!(model, :Consumption, :beta2, [beta2s[findfirst(countries .== country)] for country in dim_keys(model, :country)])

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

excel2mimi_mapping = Dict{String, Tuple{Symbol, Symbol}}("Nonlinear SAF / Distribution" => (:SAFModel, :saf_delta),
                                                         "Warming half-life / Distribution" => (:SAFModel, :FRT),
                                                         "xi_1" => (:TemperatureModel, :xi_1),
                                                         "F_2XCO2" => (:TemperatureModel, :F_2xCO2),
                                                         "T_2xCO2" => (:TemperatureModel, :fair_ECS),
                                                         "gamma" => (:TemperatureModel, :fair_gamma),
                                                         "C0" => (:TemperatureModel, :fair_C_0),
                                                         "Dataset 1" => (:CH4Model, :ch4_alpha))

function excel2mimi(model::Model, name::String; getdataframe::Function=getdataframe, tempcomp::Symbol=:TemperatureConverter)
    info = split(name, " / ")
    if name ∈ ["a_0", "a_1", "a_3", "rho_1", "rho_2", "rho_3"]
        df = getdataframe(model, :CO2Model, Symbol(replace(name, "_" => "")))
    elseif name ∈ ["r_{0) / Distribution", "r_{C} / Distribution", "r_{T} / Distribution"]
        df = getdataframe(model, :PostTemperature, Symbol(replace(info[1], "{" => "", "}" => "", ")" => "")))
    elseif name ∈ keys(excel2mimi_mapping)
        df = getdataframe(model, excel2mimi_mapping[name][1], excel2mimi_mapping[name][2])
    elseif info[2] == "Atmospheric temperature"
        df = getdataframe(model, tempcomp, :T_AT)
    elseif info[2] == "SLR total (m)"
        df = getdataframe(model, :SLRModel, :SLR)
    elseif info[2] == "World aggregate consumption per capita"
        df = getdataframe(model, :Utility, :equiv_conspc)
    elseif info[1] == "waisrate"
        df = getdataframe(model, :WAISmodel, :waisrate)
        df.waisrate *= 1000
    else
        return nothing
    end

    if length(info) == 2 && info[2] != "Distribution"
        vals = convert(Vector{Float64}, df[df.time .== parse(Int64, info[1]), 2])
    else
        vals = convert(Vector{Float64}, df[!, 2])
    end

    return vals
end
