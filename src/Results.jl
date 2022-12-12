# Load model and select configuration
include("../src/MimiMETA.jl")
model = full_model(; rcp="CP-Base")

#Update persistence parameter phi (hard-coded default: 0.5)
#myupdate_param!(model, :Consumption, :damagepersist, 0.25)

# Run the model
run(model)
explore(model) # Launches the Mimi Explorer, which is a graphical interface that lets me look at each computed component. 

#=
#= Run with MC
benchmark = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)
include("../src/lib/presets.jl")
for rr in 1:maxrr # rr number of runs
    println(rr)
    preset_fill(rr)
    run(model)
    # Save results to a data frame
end
=#

# Grab data
results = getdataframe(model, :TemperatureModel, :T_AT) # Produces initial results data frame
results[!, :SLR] = model[:SLRModel, :SLR] # Add new result variables
results[!, :total_damages_global_cumulative] = model[:TotalDamages, :total_damages_global_cumulative] # Add new result variables
#results[!, :total_damages_cumulative] = model[:TotalDamages, :total_damages_cumulative] # HOW TO ADD AN ENTIRE ARRAY RATHER THAN JUST A SINGLE CAR

#Next problem: how to grab SC-CO2 and SC-CH4 results given that they are scripts outside of the model - need to run them from here I suppose.

# Export data
using CSV
CSV.write(raw"C:\Users\Thomas\Dropbox\Tipping points in climate change economics\Mimi testbox\Results\LSE 1 Dec 2022 results.csv", results) #Default is write to working directory

#=
-SC-CH4 global
-SC-CO2 global
-SC-CH4 national
-SC-CO2 national
-total damages national

=#
=#