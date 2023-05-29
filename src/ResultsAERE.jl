##This file produces all May 2023 results for the AERE talk.

using Mimi
import Random
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")
include("../src/scch4.jl")
include("../src/scc.jl")
include("../src/bge.jl")


#=Loop over

# Scenarios
for (x,y) in [("CP-", "SSP3"), ("NP-", "SSP2"), ("1.5-", "SSP1")]]
    for z in ("Base", "GMP", "GMP-LowCH4", "GMP-HighCH4")
        rcp=x*z # Concatenate correct scenario-variant name

        ##  Load model and select configuration
        global model = full_model(;
            rcp=rcp,
            ssp=y,
    end
end

# TP configurations
for TP in ("TPs", "NoTPs")
    if TP="TPs"
        TPconfig= "interaction = true,
        pcf = "Fit of Hope and Schaefer (2016)",
        omh = "Whiteman et al. beta 20 years",
        amaz = "Cai et al. central value",
        gis = "Nordhaus central value",
        ais = "AIS",
        ism = "Value",
        amoc = "IPSL",
        nonmarketdamage = true)"
    else
        TPconfig= "interaction = false,
        pcf = "Fit of Hope and Schaefer (2016)",
        omh = "Whiteman et al. beta 20 years",
        amaz = "Cai et al. central value",
        gis = "Nordhaus central value",
        ais = "AIS",
        ism = "Value",
        amoc = "IPSL",
        nonmarketdamage = true)"

    end
end

# Phi parameters
for persistence in ("high", "default")

end
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
#=
if persistence="high"
    myupdate_param!(model, :Consumption, :damagepersist, 0.25)
end
=#

### Run the model so we can run scripts
run(model)

### Run the model in MC mode
sim_full(model, 10,
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

results = getsim_full(model, draws, true)

#=

### Calculate the balanced growth equivalent in MC mode
bgeresults = calculate_bge(model)

### Calculate the SC-CO2 in MC mode
## Miniloop over pulse year
sccresults = []
for yy in 2020:10:2100
    append!(sccresults, calculate_scc_full_mc(model,
        10, # MC reps
        "Fit of Hope and Schaefer (2016)", # PCF
        "Cai et al. central value", # AMAZ
        "Nordhaus central value", # GIS
        "none", # WAIS
        "Distribution", # SAF
        true, # ais_used
        true, # ism_used
        true, # omh_used
        true, # amoc_used
        false, # persist
        false, # emuc
        false, # prtp
        yy, # pulse year
        10.0, # pulse size
        1.5)) # EMUC
end

### Calculate the SC-CH4 in MC mode
## Miniloop over pulse year
scch4results = []
for yy in 2020:10:2100
    append!(scch4results, calculate_scch4_full_mc(model,
        10, # MC reps
        "Fit of Hope and Schaefer (2016)", # PCF
        "Cai et al. central value", # AMAZ
        "Nordhaus central value", # GIS
        "none", # WAIS
        "Distribution", # SAF
        true, # ais_used
        true, # ism_used
        true, # omh_used
        true, # amoc_used
        false, # persist
        false, # emuc
        false, # prtp
        yy, # pulse year
        0.36, # pulse size
        1.5)) # EMUC
end



##Settings
#-phi=0.5 and phi=0.25
#-SAF on always
#-Non-market damages on always



##Results post Julia
#-Truncate runs as in PNAS paper
#-Write Stata script to make pretty graphs and maps


#=Monte Carlo to do post AERE
#-Ensure individual runs are comparable within MC draw
=#

=#