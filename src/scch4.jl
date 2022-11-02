using Mimi

function calculate_scch4(model::Model, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    mm = create_marginal_model(model, pulse_size)

    pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)
    run(mm)

    mm.modified[:CH4Model, :ch4_extra][pulse_index] = pulse_size
    run(mm)

    globalwelfare_marginal = sum(mm[:Utility, :world_disc_utility][pulse_index:191])

    global_conspc = sum(mm.base[:Consumption, :conspc][pulse_index, :] .* mm.base[:Utility, :pop][pulse_index, :]) / mm.base[:Utility, :world_population][pulse_index]
    -(globalwelfare_marginal / (global_conspc^-emuc)) / 1e6 #Copying META scaling since emissions of CO2 are in Gt but CH4 in Mt.
end

include("../src/MimiMETA.jl")
model = base_model(; tdamage="pointestimate", slrdamage="mode")

#Try different pulse sizes starting from the upper bound of the equivalent to META's 10GtCO2 pulse
calculate_scch4(model, 2020, 0.36, 1.5)#360,000 tCH4
calculate_scch4(model, 2020, 0.24, 1.5)#240,000 tCH4
calculate_scch4(model, 2020, 0.12, 1.5)#120,000 tCH4
#calculate_scch4(model, 2020, 0.06, 1.5)#60,000 tCH4
#calculate_scch4(model, 2020, 0.03, 1.5)#30,000 tCH4
#calculate_scch4(model, 2020, 0.015, 1.5)#15,000 tCH4
#calculate_scch4(model, 2020, 0.01, 1.5)#10,000 tCH4
#calculate_scch4(model, 2020, 0.005, 1.5)#5,000 tCH4
#calculate_scch4(model, 2020, 0.0025, 1.5)#2,500 tCH4
#calculate_scch4(model, 2020, 0.002, 1.5)#2,000 tCH4
#calculate_scch4(model, 2020, 0.001, 1.5)#1,000 tCH4
#calculate_scch4(model, 2020, 0.0005, 1.5)#500 tCH4
#calculate_scch4(model, 2020, 0.00005, 1.5) #50 tCH4
#calculate_scch4(model, 2020, 0.000025, 1.5) #25 tCH4
#calculate_scch4(model, 2020, 0.000005, 1.5) #5 tCH4
#calculate_scch4(model, 2020, 0.000001, 1.5) #1 tCH4

#calculate_scch4(model, 2020, 2500., 1.5)
#calculate_scch4(model, 2020, 2000., 1.5)
#calculate_scch4(model, 2020, 1000., 1.5)
#calculate_scch4(model, 2020, 500., 1.5)
#calculate_scch4(model, 2020, 0.000001, 1.5) #1 tCH4
#Result:
#-SC-CH4 extremely stable: 1151.36USD/tCH4 for a 1 tCH4 pulse, and 1150.57 for a 360,000 tCH4 pulse. Need to check with TPs later.
#Problem:
#-VisualStudio's terminal only shows the output of the last line of code. 
#-VisualStudio's terminal is slow
