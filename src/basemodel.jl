using Mimi
using DataFrames, CSV
countries = CSV.read("../data/pattern-scaling_new.csv", DataFrame).Country
include("lib/gdppc.jl")

function test_model(; startyear::Int64=2010)
    m = Model()

    set_dimension!(m, :time, collect(startyear:2200))
    set_special_model_dimensions!(m)

    return m
end

function set_special_model_dimensions!(m::Model)
    set_dimension!(m, :country, countries)
    set_dimension!(m, :region, unique(gdps_ssp.REGION_SHORT))
    set_dimension!(m, :monsoonsteps, 35)
end

