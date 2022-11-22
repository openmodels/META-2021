using Mimi
using DataFrames, CSV
countries = CSV.read("../data/pattern-scaling.csv", DataFrame).Country
include("lib/gdppc.jl")

function test_model()
    m = Model()

    set_dimension!(m, :time, collect(2020:2200)) #I start the model only in 2020 in line with our new scenario data.
    set_dimension!(m, :country, countries)
    set_dimension!(m, :region, unique(gdps_ssp.REGION_SHORT))
    set_dimension!(m, :monsoonsteps, 35)

    return m
end
