using Mimi

function calculate_scc(model::Model, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    mm = create_marginal_model(model, pulse_size)

    pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)
    run(mm)

    mm.modified[:CO2Model, :co2_extra][pulse_index] = pulse_size
    run(mm)

    globalwelfare_marginal = sum(mm[:Utility, :world_disc_utility][pulse_index:191])

    global_conspc = sum(mm.base[:Consumption, :conspc][pulse_index, :] .* mm.base[:Utility, :pop][pulse_index, :]) / mm.base[:Utility, :world_population][pulse_index]
    -(globalwelfare_marginal / (global_conspc^-emuc)) / 1e9
end

include("../src/MimiMETA.jl")
model = base_model(; tdamage="pointestimate", slrdamage="mode")
calculate_scc(model, 2020, 10., 1.5)
