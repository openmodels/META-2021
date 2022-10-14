using Test, CSV, DataFrames
include("../src/MimiMETA.jl")
import Mimi.has_parameter

benchmark = CSV.read("../data/benchmark/ExcelMETA-notp.csv", DataFrame)

## Setup the model

function myupdate_param!(model, comp, param, value)
    unique_name = Symbol("$(comp)_$param")
    if has_parameter(model.md, unique_name)
        update_param!(model, unique_name, value)
    else
        set_param!(model, comp, param, unique_name, value)
    end
end

for rr in 1:nrow(benchmark)
    global model = base_model() # XXX: Consumption currently refers to this...

    ## Fill in values
    beta1indexes = findall(x -> occursin("beta1dist", x), names(benchmark))
    countries = [x[1:3] for x in names(benchmark)[beta1indexes]]
    beta1s = convert(Array, benchmark[rr, beta1indexes])
    beta2indexes = findall(x -> occursin("beta2dist", x), names(benchmark))
    beta2s = convert(Array, benchmark[rr, beta2indexes])

    myupdate_param!(model, :Consumption, :seeds, zeros(dim_count(model, :country)))
    myupdate_param!(model, :Consumption, :beta1, [beta1s[findfirst(countries .== country)] for country in dim_keys(model, :country)])
    myupdate_param!(model, :Consumption, :beta2, [beta2s[findfirst(countries .== country)] for country in dim_keys(model, :country)])
    
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
    myupdate_param!(model, :Forcing, :F_2xCO2, benchmark.F_2XCO2[rr])
    myupdate_param!(model, :CH4Model, :ch4_alpha, benchmark."Dataset 1"[rr])

    myupdate_param!(model, :Consumption, :damagepersist, 0.5)

    slrindexes = findfirst(names(benchmark) .== "Distribution / AFG"):findfirst(names(benchmark) .== "Distribution / ZWE")
    countries = [x[end-2:end] for x in names(benchmark)[slrindexes]]
    slrcoeffs = convert(Array, benchmark[rr, slrindexes])
    myupdate_param!(model, :Consumption, :slrcoeff, [slrcoeffs[findfirst(countries .== country)] for country in dim_keys(model, :country)])

    myupdate_param!(model, :Utility, :EMUC, 1.5)
    myupdate_param!(model, :Utility, :PRTP, 0.01)

    run(model)

    ## Test the model

    T_AT = model[:TemperatureModel, :T_AT][11:10:191]
    T_AT_compare = convert(Array, benchmark[rr, 2:20])

    @test T_AT ≈ T_AT_compare

    SLR = model[:SLRModel, :SLR][11:10:191]
    SLR_compare = convert(Array, benchmark[rr, 21:39])

    @test SLR ≈ SLR_compare

    globalwelfare = sum(model[:Utility, :world_disc_utility][11:191])
    globalwelfare_compare = benchmark."Global welfare"[rr]

    @test globalwelfare ≈ globalwelfare_compare
end
