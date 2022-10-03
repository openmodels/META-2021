using Test, CSV
include("../src/basemodel.jl")

benchmark = CSV.read("../data/benchmark/GISModel.csv", DataFrame)

## Setup the model

model = test_model()
include("../src/components/GIS.jl")

gismodel = addGISModel(model, "Nordhaus central value")
gismodel[:T_AT] = benchmark."Atmospheric temperature" 

run(model)

## Test the model

SLR_GIS = model[:GISModel, :SLR_GIS]
SLR_GIS_compare = benchmark."SLR_GIS"

@test SLR_GIS ≈ SLR_GIS_compare rtol=1e-4
