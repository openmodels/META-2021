using CSV, Tables
include("../src/MimiMETA.jl")
include("../src/lib/presets.jl")

## BASE RUNS WITHOUT EXTRA TP: ONLY WAIS, GIS, SAF
for (x,y) in [("CP-", "SSP3"), ("NP-", "SSP2"), ("1.5-", "SSP1")]#, ("1.5-", "SSP2")]
    for z in ("Base", "GMP", "GMP-LowCH4", "GMP-HighCH4") 
        rcp=x*z # Concatenate correct scenario-variant name
              
        ##  Load model and select configuration
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
            pcf="Fit of Hope and Schaefer (2016)",
            omh="Whiteman et al. beta 20 years",
            amaz="Cai et al. central value",
            gis="Nordhaus central value", 
            wais="Value", 
            ism="Value", 
            amoc="IPSL", 
            nonmarketdamage=false)

        ## Update persistence parameter phi (hard-coded default: 0.5)
        myupdate_param!(model, :Consumption, :damagepersist, 0.25)

        ## Run the model
        run(model)
        #explore(model) # Launches the Mimi Explorer, which is a graphical interface that lets me look at each computed component. 
        subres = []
        for yy in 2020:10:2100
            ## Post-model scripts for SC-X for pulses in 2020
            include("../src/scch4.jl")
            append!(subres, calculate_scch4(model, yy, 0.001, 1.5))
            #sc_ch4_mc = calculate_scch4_mc(model, preset_fill, nrow(benchmark), 2020, 0.00001, 1.5) NEED TO ADD CODE TO INIT MC RUNS

            include("../src/scc.jl")
            append!(subres, calculate_scc(model, yy, 10., 1.5))
            #sc_co2_mc = calculate_scc_mc(model, preset_fill, nrow(benchmark), 2020, 10., 1.5) NEED TO ADD CODE TO INIT MC RUNS

            append!(subres, yy)
        end

        ## Append to full results data
        append!(subres, rcp)

        global results = [results; subres]
               
    end
end

## Export data
CSV.write("SC-X_allTPs.csv", Tables.table(results), writeheader=false)

## Clean results variables
results=nothing
