using Mimi
using DataFrames, CSV

include("../src/MimiMETA.jl")

model = base_model(; rcp="CP-GMP", tdamage="pointestimate", slrdamage="mode")
run(model)

function calculate_bge(model::Model, outpath::String="")
    results = DataFrame(country=String[], bge=Float64[])

    emuc = model[:BGE, :EMUC]

    # Present value of discounted utility streams
    sum_world_disc_utility                  = sum(model[:BGE, :world_disc_utility][1:191])
    sum_world_disc_utility_counterfactual   = sum(model[:BGE, :world_disc_utility_counterfactual][1:191])

    # Prefill one dimensional array of zeros for each country
    isos = dim_keys(model, :country)
    sum_disc_utility                        = zeros(length(isos), 1)
    sum_disc_utility_counterfactual         = zeros(length(isos), 1)

    for cc in 1:length(isos)
        sum_disc_utility[cc]                = sum(model[:BGE, :disc_utility][1:191, cc])
        sum_disc_utility_counterfactual[cc] = sum(model[:BGE, :disc_utility_counterfactual][1:191, cc])
    end

    # Apply BGE method (eq. 5 in Anthoff and Tol 2009 ERE)
    world_bge = (sum_world_disc_utility_counterfactual/sum_world_disc_utility)^(1/(1-emuc))-1
    push!(results, ("globe", world_bge))

    bge                                     = zeros(length(isos), 1)
    for cc in 1:length(isos)
        bge[cc] = (sum_disc_utility_counterfactual[cc]/sum_disc_utility[cc])^(1/(1-emuc))-1
        push!(results, (dim_keys(model, :country)[cc], bge[cc]))
    end

    if outpath != ""
        CSV.write(outpath, results)
    end

    results
end

calculate_bge(model, outpath="bges-CP-GMP-1.5.csv")

results = DataFrame(country=String[], bge=Float64[], rcp=String[])

for rcp in ["NP-Base", "CP-Base", "1.5-Base", "NP-GMP", "CP-GMP", "1.5-GMP", "NP-GMP-LowCH4", "CP-GMP-LowCH4", "1.5-GMP-LowCH4", "NP-GMP-HighCH4", "CP-GMP-HighCH4", "1.5-GMP-HighCH4"]
    model = base_model(; rcp=rcp, tdamage="pointestimate", slrdamage="mode")
    run(model)
    subres = calculate_bge(model)
    subres[!, :rcp] .= rcp

    results = [results; subres]
end

CSV.write("bges.csv", results)
