##This file produces all May 2023 results for the AERE talk.

using Mimi
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl") # Do we need presets.jl or does montecarlo.jl supersede that?
include("../src/scch4.jl")
include("../src/scc.jl")
include("../src/bge.jl")

## Set common MC seed (might need to redo the same for scch4.jl and scc.jl)

## Set number of MC repetitions
reps = 10 #Use 10000 for actual results once code is stable

#=QUESTIONS
-montecarlo.jl already sets a 'trials' variable that is the number of reps - how to access that here?
-More generally, what is the connection between MimiMETA.jl, and montecarlo.jl?
=#

## Settings for model
mm = (model, reps,
    "Fit of Hope and Schaefer (2016)", # PCF
    "Cai et al. central value", # AMAZ
    "Nordhaus central value", # GIS
    "none", # WAIS
    "Distribution", # SAF
    true, # AIS
    true, # ISM
    true, # OMH
    true, # AMOC
    false, # persist
    false, # emuc
    false # prtp
)

#=QUESTIONS
-Can't find an option to set non-market damages in montecarlo.jl anymore
-Where is EMUC's default set?
-Where is PRTP's default set? 
=#

## Compile model
run(mm)

##Settings
#-phi=0.5
#-SAF on always
#-Non-market damages on always

##Loop over
#-SAF only, all TPs (AIS for Antarctica) on
#-Scenarios: 3 base and 3 RMA scenarios




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
