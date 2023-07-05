using Test
include("../src/MimiMETA.jl")

model_base = base_model(rcp="RCP4.5", ssp="SSP2", tdamage="pointestimate", slrdamage="mode")
run(model_base)

conspc_base = model_base[:Consumption, :conspc][end, dim_keys(model_base, :country) .== "GBR"]

model_base2 = full_model(rcp="RCP4.5", ssp="SSP2", saf=false, interaction=false, pcf=false, omh=false, amaz=false, gis=false, ais=false, ism=false, amoc=false)
run(model_base2)

conspc_base2 = model_base2[:Consumption, :conspc][end, dim_keys(model_base2, :country) .== "GBR"]

@test conspc_base == conspc_base2

model_amoc = full_model(rcp="RCP4.5", ssp="SSP2", saf=false, interaction=false, pcf=false, omh=false, amaz=false, gis=false, ais=false, ism=false, amoc="IPSL")
run(model_amoc)

conspc_amoc = model_amoc[:Consumption, :conspc][end, dim_keys(model_amoc, :country) .== "GBR"]

@test conspc_base != conspc_amoc

