using Mimi
using DataFrames, CSV

include("../src/MimiMETA.jl")
include("../src/lib/presets.jl")

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

#calculate_bge(model, "bges-CP-GMP-1.5.csv") # Not needed, but keeping this in case it was a placeholder for James for later

global results = DataFrame(country=String[], bge=Float64[], rcp=String[])

for (x,y) in [("CP-", "SSP2"), ("NP-", "SSP3"), ("1.5-", "SSP1")]#, ("1.5-", "SSP2")] # Loop over base scenarios and socioeconomics
    for z in ("Base", "GMP", "GMP-LowCH4", "GMP-HighCH4") # Loop over CH4 variants
        rcp=x*z # Concatenate correct scenario-variant name
            
        ##  Load model and select configuration
        include("../src/MimiMETA.jl")
        model = full_model(; 
            rcp=rcp, 
            ssp=y, 
            co2="Expectation", 
            ch4="default", 
            warming="Best fit multi-model mean", 
            tdamage="pointestimate", 
            slrdamage="mode", 
            saf="Distribution mean", 
            interaction=true, 
            pcf=false,#"Fit of Hope and Schaefer (2016)", 
            omh=false,#"Whiteman et al. beta 20 years", 
            amaz=false,#"Cai et al. central value", 
            gis="Nordhaus central value", 
            wais="Value", 
            ism=false,#"Value", 
            amoc=false,#"IPSL", 
            nonmarketdamage=false) 
        
        # Set persistence parameter to 0.25
        myupdate_param!(model, :Consumption, :damagepersist, 0.25)

        #Other robustness
        #myupdate_param!(model, :BGE, :PRTP, 0.02) 
        
        run(model)
        subres = calculate_bge(model)
        subres[!, :rcp] .= rcp

        global results = [results; subres]
    end
end

CSV.write("bges.csv", results)
