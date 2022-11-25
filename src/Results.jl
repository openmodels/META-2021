# Load model
using Mimi
using MimiMETA
m = MimiMETA.get_model()

///Configure the model - although I believe this is part of the step that loads the model

*Change parameter values
set_param!(m, :tcr_transientresponse, 3) #This is still MimiPAGE readme code.

///Run the model
run(m)

///Extract and store variables of interest
result.m #returns the Model: an array full of different variables. How can I view the Mimi results viewer to decide what I need? (workaround: just look at MimiMETA component .jl files directly)

marginal_temp = result.mm[:ClimateTemperature, :rt_realizedtemperature] #MimiPAGE code: Mimi seems to organize the array by modules/sub-arrays: find module name, extract the variable of interest.