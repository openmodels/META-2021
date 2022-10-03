using Test, CSV, Mimi, Random, Distributions
include("../src/basemodel.jl")
include("../src/components/WAIS.jl")

benchmark = CSV.read("../data/benchmark/WAIS.csv", DataFrame) #Generated with RCP8.5

## Set up the model

model = test_model()
WAISmodel = addWAISmodel(model, "Value")

#Add temperature time series from benchmark file
WAISmodel[:T_AT] = benchmark."Temperature" 
#For the non-Monte Carlo test, only compare cumulative WAIS SLR given a sequence of WAIS indicator function realizations.
WAISmodel[:uniforms] = 1 .-benchmark."Binomial draw"
WAISmodel[:waisrate] = 0.027112639

run(model)

## Test the model

SLR_WAIS = model[:WAISmodel, :SLR_WAIS]
SLR_WAIS_compare = benchmark."Cumulative SLR_WAIS (m)"

@test SLR_WAIS ≈ SLR_WAIS_compare rtol=1e-4