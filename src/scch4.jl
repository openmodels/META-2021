using Mimi

function calculate_scch4(model::Model, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    mm = create_marginal_model(model, pulse_size)

    pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)
    run(mm)

    mm.modified[:CH4Model, :ch4_extra][pulse_index] = pulse_size
    run(mm)

    globalwelfare_marginal = sum(mm[:Utility, :world_disc_utility][pulse_index:191])

    global_conspc = sum(mm.base[:Consumption, :conspc][pulse_index, :] .* mm.base[:Utility, :pop][pulse_index, :]) / mm.base[:Utility, :world_population][pulse_index]
    -(globalwelfare_marginal / (global_conspc^-emuc)) / 1e6 
end

include("../src/MimiMETA.jl")
model = base_model(; rcp="NP-Base", tdamage="pointestimate", slrdamage="mode")

calculate_scch4(model, 2020, 0.36, 1.5)#360,000 tCH4 pulse

function calculate_scch4_mc(model::Model, preset_fill::Function, maxrr::Int64, pulse_year::Int64, pulse_size::Float64, emuc::Float64)
    scch4s = []
    for rr in 1:maxrr
        println(rr)
        preset_fill(rr)
        push!(scch4s, calculate_scch4(model, pulse_year, pulse_size, emuc))
    end
    scch4s
end

include("../src/lib/presets.jl")
benchmark = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)
model = full_model()
preset_fill(rr) = preset_fill_tp(model, benchmark, rr)
calculate_scch4_mc(model, preset_fill, 5, 2020, 0.36, 1.5) # Runs 500 MC reps.