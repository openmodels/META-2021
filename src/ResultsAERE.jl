##This file produces all May 2023 results for the AERE talk.

using Mimi
import Random
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")
## include("../src/lib/presets.jl") # Needed?
include("../src/scch4.jl")
include("../src/scc.jl")
include("../src/bge.jl")


#=Loop over

# Scenarios
6 SCENARIOS

# TP configurations
2 TP CONFIGS

# Phi parameters
2 PHI SETTINGS
=#

Random.seed!(26052023)

### Create the model (this just sets a long string that the other functions can use)
model = full_model(;
    rcp = "CP-Base",
    ssp = "SSP2",
    co2 = "Expectation",
    ch4 = "default",
    warming = "Best fit multi-model mean",
    tdamage = "pointestimate",
    slrdamage = "mode",
    saf = "Distribution mean",
    interaction = true,
    pcf = "Fit of Hope and Schaefer (2016)",
    omh = "Whiteman et al. beta 20 years",
    amaz = "Cai et al. central value",
    gis = "Nordhaus central value",
    ais = "AIS",
    ism = "Value",
    amoc = "IPSL",
    nonmarketdamage = true)

### Update the persistence parameter, or not


### Run the model in MC mode
sim_full(model, 100,
         "Fit of Hope and Schaefer (2016)", # PCF
         "Cai et al. central value", # AMAZ
         "Nordhaus central value", # GIS
         "none", # WAIS
         "Distribution", # SAF
         true, # ais_used
         true, # ism_used
         true, # omh_used
         true, # amoc_used
         false, # persit
         false, # emuc
         false) # prtp
#getsim_full(inst::Union{ModelInstance, MarginalInstance}, draws::DataFrame; save_rvs::Bool=true)

### Calculate the balanced growth equivalent in MC mode
run(model)
bgeresults = calculate_bge(model)

### Calculate the SC-CO2 in MC mode
sccresults = calculate_scc_full_mc(model, 100, 10.0, 1.5)
sccresults = calculate_scc(model, 2020, 10.0, 1.5)
#Miniloop over 2020(10)2100

### Calculate the SC-CH4 in MC mode
scch4results = calculate_scch4(model, 2020, 0.36, 1.5)
#Miniloop over 2020(10)2100



#=Outstanding to do
-Set common MC seed (might need to redo the same for scch4.jl and scc.jl)=#





##Settings
#-phi=0.5 and phi=0.25
#-SAF on always
#-Non-market damages on always




scch4s = calculate_scch4_full_mc(model, 100,
                              "Fit of Hope and Schaefer (2016)", # PCF
                              "Cai et al. central value", # AMAZ
                              "Nordhaus central value", # GIS
                              "Distribution", # WAIS
                              "Distribution", # SAF
                              false, # ais_used
                              true, # ism_used
                              false, # omh_used
                              true, # amoc_used
                              false, # persit
                              false, # emuc
                              false, # prtp
                              2020, 0.36, 1.5)


##Results post Julia
#-Truncate runs as in PNAS paper
#-Write Stata script to make pretty graphs and maps


#=Monte Carlo to do post AERE
#-Ensure individual runs are comparable within MC draw
=#
