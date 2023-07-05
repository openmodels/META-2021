using CSV
using Mimi
include("../src/MimiMETA.jl")

model = base_model(; rcp="CP-Base", tdamage="pointestimate", slrdamage="mode")
run(model)
df = getdataframe(model, :temperature, :T)

CSV.write("orig.csv", df)

model = base_model(; rcp="CP-Base", tdamage="pointestimate", slrdamage="mode", useexrf=true)
run(model)
df = getdataframe(model, :temperature, :T)

CSV.write("useex.csv", df)
