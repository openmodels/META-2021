using Test, CSV, DataFrames
include("../src/MimiMETA.jl")
include("../src/scc.jl")
include("../src/lib/presets.jl")

benchmark_notp = CSV.read("../data/benchmark/ExcelMETA-notp.csv", DataFrame)
benchmark_tp = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)

rcps = ["RCP3-PD/2.6", "RCP4.5", "RCP6", "RCP8.5", "RCP4.5"]
ssps = ["SSP1", "SSP2", "SSP4", "SSP5", "SSP5"]
scc_notp = [33.96, 52.03, 80.55, 32.85, 23.12]
scc_tp = [45.42, 64.80, 92.91, 39.23, 28.90]

mine_scc_notp = []
mine_scc_tp = []
for ii in 1:length(rcps)
    println("$(rcps[ii]), $(ssps[ii])")
    model_notp = base_model(; rcp=rcps[ii], ssp=ssps[ii])
    model_tp = model = full_model(; rcp=rcps[ii], ssp=ssps[ii])

    global model = model_notp
    sccs_notp = calculate_scc_mc(model, (rr) -> preset_fill_notp(model, benchmark_notp, rr), nrow(benchmark_notp), 2020, 10., 1.5)

    # @test mean(sccs_notp) ≈ scc_notp[ii] atol=1e-1
    push!(mine_scc_notp, mean(sccs_notp))

    global model = model_tp
    sccs_tp = calculate_scc_mc(model, (rr) -> preset_fill_tp(model, benchmark_tp, rr), nrow(benchmark_tp), 2020, 10., 1.5)

    # @test mean(sccs_tp) ≈ scc_tp[ii] atol=1e-1
    push!(mine_scc_tp, mean(sccs_tp))
end

(mine_scc_tp .- mine_scc_notp) ./ mine_scc_notp
