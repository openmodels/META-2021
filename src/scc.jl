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

model = base_model(; rcp="CP-Base", tdamage="pointestimate", slrdamage="mode")
#calculate_scc(model, 2020, 10., 1.5)

function calculate_scc_mc(model::Model, preset_fill::Function, maxrr::Int64, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    sccs = []
    for rr in 1:maxrr
        println(rr)
        preset_fill(rr)
        push!(sccs, calculate_scc(model, pulse_year, pulse_size, emuc))
    end
    sccs
end

include("../src/lib/presets.jl")
benchmark = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)
model = full_model()
preset_fill(rr) = preset_fill_tp(model, benchmark, rr)
#calculate_scc_mc(model, preset_fill, nrow(benchmark), 2020, 10., 1.5) # Runs 500 MC reps.
