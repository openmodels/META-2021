using Test, CSV, Mimi, Random, Distributions, XLSX

include("../src/basemodel.jl")
include("../src/components/AIS.jl")

benchmark = CSV.read("../data/benchmark/AIS.csv", DataFrame)

## Set up the model

model = test_model()
AISmodel = addAISmodel(model, "Value")

#Add temperature time series from benchmark file

temperature_data = XLSX.readdata("../data/Meta model for AIS melting May 2022.xlsx", "EAIS Basal Melt", "B112:B302")
temperature_data_1900 = XLSX.readdata("../data/Meta model for AIS melting May 2022.xlsx", "EAIS Basal Melt", "B2:B302")

Temps = CSV.read("../data/Temperature since 1900.csv", DataFrame)
R_functions_EAIS = XLSX.readdata("../data/Response functions.xlsx", "EAIS", "B3:R193")
R_functions_Ross = XLSX.readdata("../data/Response functions.xlsx", "Ross", "B3:R193") 
R_functions_Amundsen = XLSX.readdata("../data/Response functions.xlsx", "Amundsen", "B3:R193") 
R_functions_Weddell = XLSX.readdata("../data/Response functions.xlsx", "Weddell", "B3:R193") 
R_functions_Peninsula = XLSX.readdata("../data/Response functions.xlsx", "Peninsula", "B3:R193")
basal_melt_EAIS = XLSX.readdata("../data/basal_melt_models.xlsx", "EAIS", "A2:C20")
basal_melt_Ross = XLSX.readdata("../data/basal_melt_models.xlsx", "Ross", "A2:C20")
basal_melt_Amundsen = XLSX.readdata("../data/basal_melt_models.xlsx", "Amundsen", "A2:C20")
basal_melt_Weddell = XLSX.readdata("../data/basal_melt_models.xlsx", "Weddell", "A2:C20")
basal_melt_Peninsula = XLSX.readdata("../data/basal_melt_models.xlsx", "Peninsula", "A2:C20")

AISmodel[:T_AT] = Temps[111:301, "Global Mean Surface Temperature"]
AISmodel[:temp_1900] = Temps[:, "Global Mean Surface Temperature"]
AISmodel[:R_functions_EAIS] = R_functions_EAIS[:, 9]
AISmodel[:R_functions_Ross] = R_functions_Ross[:, 9]
AISmodel[:R_functions_Amundsen] = R_functions_Amundsen[:, 9]
AISmodel[:R_functions_Weddell] = R_functions_Weddell[:, 9]
AISmodel[:R_functions_Peninsula] = R_functions_Peninsula[:, 9]
AISmodel[:β_EAIS] = basal_melt_EAIS[10, 2]
AISmodel[:δ_EAIS] = basal_melt_EAIS[10, 3]
AISmodel[:β_Ross] = basal_melt_Ross[10, 2]
AISmodel[:δ_Ross] = basal_melt_Ross[10, 3]
AISmodel[:β_Amundsen] = basal_melt_Amundsen[10, 2]
AISmodel[:δ_Amundsen] = basal_melt_Amundsen[10, 3]
AISmodel[:β_Weddell] = basal_melt_Weddell[10, 2]
AISmodel[:δ_Weddell] = basal_melt_Weddell[10, 3]
AISmodel[:β_Peninsula] = basal_melt_Peninsula[10, 2]
AISmodel[:δ_Peninsula] = basal_melt_Peninsula[10, 3]

run(model)

## Test the model

SLR_AIS = model[:AISmodel, :SLR_AIS]
SLR_AIS_compare = benchmark."Cumulative Dynamic and SMB SLR"

@test SLR_AIS ≈ SLR_AIS_compare rtol=1e-4