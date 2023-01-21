using Test, CSV
include("../src/basemodel.jl")
include("../src/components/RCP.jl")

## Setup the model

model = test_model()
template = addRCP(model, "RCP8.5")

run(model)
