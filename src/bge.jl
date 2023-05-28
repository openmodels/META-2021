using Mimi
using DataFrames, CSV

include("../src/MimiMETA.jl")
include("../src/lib/presets.jl")

#model = base_model(; rcp="CP-GMP", tdamage="pointestimate", slrdamage="mode")
#run(model)

function calculate_bge(model::Model, outpath::String="")
    results = DataFrame(country=String[], bge=Float64[])

    emuc = model[:BGE, :EMUC]
    startyear = findfirst(dim_keys(model, :time) .== 2010)

    # Present value of discounted utility streams
    sum_world_disc_utility                  = sum(model[:BGE, :world_disc_utility][startyear:end])
    sum_world_disc_utility_counterfactual   = sum(model[:BGE, :world_disc_utility_counterfactual][startyear:end])

    # Prefill one dimensional array of zeros for each country
    isos = dim_keys(model, :country)
    sum_disc_utility                        = zeros(length(isos), 1)
    sum_disc_utility_counterfactual         = zeros(length(isos), 1)

    for cc in 1:length(isos)
        sum_disc_utility[cc]                = sum(model[:BGE, :disc_utility][startyear:end, cc])
        sum_disc_utility_counterfactual[cc] = sum(model[:BGE, :disc_utility_counterfactual][startyear:end, cc])
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
