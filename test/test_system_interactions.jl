using Test, CSV, DataFrames
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")

allres = DataFrame(aiscalib=String[], useinteract=Bool[], slr2200=Float64[], stderr=Float64[])
for aiscalib in ["WAIS", "AIS"]
    for useinteract in [false, true]
        global model = full_model(rcp="RCP4.5", ssp="SSP2", ais=aiscalib, interaction=useinteract)
        run(model)

        results = sim_full(model, 100, "Fit of Hope and Schaefer (2016)", # PCF
                           "Cai et al. central value", # AMAZ
                           "Nordhaus central value", # GIS
                           (aiscalib == "WAIS" ? "Distribution" : "none"), # WAIS
                           "Distribution", # SAF
                           (aiscalib == "AIS"), # ais_used
                           true, # ism_used
                           true, # omh_used
                           true, # amoc_used
                           false, # persit
                           false, # emuc
                           false;
                           getsim=(inst, draws; save_rvs) -> inst[:SLRModel, :SLR][end]) # prtp
        push!(allres, [aiscalib, useinteract, mean(results[:other]), std(results[:other]) / sqrt(length(results[:other]))])
    end
end
