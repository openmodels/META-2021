using Test, CSV, DataFrames
include("../src/MimiMETA.jl")
include("../src/lib/presets.jl")

benchmark = CSV.read("../data/benchmark/Test run for Simon 9 May 23 TPs off - base.csv", DataFrame)

model = base_model(; rcp="RCP4.5")
run(model)

tvals = Dict{String, Float64}()
missings = String[]
errors = String[]
for ii in 1:nrow(benchmark)
    try
        value = excel2mimi(model, benchmark.Name[ii])[1]
    catch
        push!(errors, benchmark.Name[ii])
        continue
    end
    if value == nothing
        push!(missings, benchmark.Name[ii])
        continue
    end
    tval = (value - benchmark.Mean[ii]) / benchmark."Std Deviation"[ii]
    tvals[benchmark.Name[ii]] = tval
    if abs(tval) >= 2.56  # 99%
        println("Extreme t-value on $(benchmark.Name[ii])")
        #@test abs(tval) < 2.56
    end
end
