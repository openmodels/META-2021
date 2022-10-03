using Test, CSV
include("../src/basemodel.jl")
include("../src/components/Template.jl")

## Setup the model

model = test_model()
template = addTemplate(model)

template[:myinput] = repeat([1], dim_count(model, :time))

run(model)

## Test the model

myvar = model[:Template,  :myvar]
myvar_compare = collect(101:dim_count(model, :time)+100)

@test myvar ≈ myvar_compare rtol=1e-4
