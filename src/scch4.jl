using Mimi

function calculate_scch4(model::Model, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    mm = create_marginal_model(model, pulse_size)

    pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)
    run(mm)

    mm.modified[:CH4Model, :ch4_extra][pulse_index] = pulse_size
    run(mm)

    globalwelfare_marginal = sum(mm[:Utility, :world_disc_utility][pulse_index:181])

    global_conspc = sum(mm.base[:Consumption, :conspc][pulse_index, :] .* mm.base[:Utility, :pop][pulse_index, :]) / mm.base[:Utility, :world_population][pulse_index]
    -(globalwelfare_marginal / (global_conspc^-emuc)) / 1e6 
end

include("../src/MimiMETA.jl")
model = base_model(; rcp="CP-Base", tdamage="pointestimate", slrdamage="mode")

calculate_scch4(model, 2021, 0.36, 1.5)#360,000 tCH4 pulse
#Once uncertainty is fully implemented: try different pulse sizes starting from the upper bound of the equivalent to META's 10GtCO2 pulse, 360,000 tCH4.
