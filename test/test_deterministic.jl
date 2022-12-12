include("../src/MimiMETA.jl")

for ii in 1:5
    global model = full_model(;
                       rcp="RCP3-PD/2.6",
                       ssp="SSP1",
                       co2="Expectation",
                       ch4="default",
                       warming="Best fit multi-model mean",
                       tdamage="pointestimate",
                       slrdamage="mode",
                       saf=false,#"Distribution mean",
                       interaction=true,
                       pcf=false,#"Fit of Hope and Schaefer (2016)",
                       omh=false,#"Whiteman et al. beta 20 years",
                       amaz=false,#"Cai et al. central value",
                       gis=false,#"Nordhaus central value",
                       wais=false,#"Value",
                       ism=false,#"Value",
                       amoc=false,#"IPSL",
                       nonmarketdamage=false)

    run(model)

    include("../src/scc.jl")
    sc_co2_global = calculate_scc(model, 2020, 10., 1.5)

    println(sc_co2_global)
end
