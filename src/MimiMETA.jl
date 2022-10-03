using Distributions

include("../src/basemodel.jl")
include("../src/components/RCP.jl")
include("../src/components/CO2Model.jl")
include("../src/components/CH4Model.jl")
include("../src/components/Forcing.jl")
include("../src/components/SAF.jl")
include("../src/components/TemperatureModel.jl")
include("../src/components/PostTemperature.jl")
include("../src/components/PCF.jl")
include("../src/components/OMH.jl")
include("../src/components/AmazonDieback.jl")
include("../src/components/GIS.jl")
include("../src/components/WAIS.jl")
include("../src/components/SLR.jl")
include("../src/components/ISM.jl")
include("../src/components/PatternScaling.jl")
include("../src/components/AMOC.jl")
include("../src/components/Interactions.jl")
include("../src/components/Consumption.jl")
include("../src/components/NonMarketDamages.jl")
include("../src/components/Utility.jl")

function base_model()
    model = test_model();

    RCPmodel = addRCP(model, "RCP8.5");
    co2model = addCO2Model(model, "Expectation");
    CH4model = addCH4Model(model, "Value");
    forcing = addForcing(model, "Best fit multi-model mean");
    temperaturemodel = addTemperatureModel(model, "Best fit multi-model mean");
    posttemp = addPostTemperature(model, "Expectation");
    slr = addSLR(model);
    pattscale = addPatternScaling(model);
    cons = addConsumption(model, "none", "none", "SSP2");
    utility = addUtility(model, "SSP2");

    # Setup CO2 model
    co2model[:co2_rcp] = RCPmodel[:co2_rcp];
    co2model[:alpha] = posttemp[:alpha];

    # Setup CH4 model
    CH4model[:ch4_rcp] = RCPmodel[:ch4_rcp];
    CH4model[:ch4_conc_rcp] = RCPmodel[:ch4_conc_rcp];
    CH4model[:n2o_conc_rcp] = RCPmodel[:n2o_conc_rcp];

    # Setup forcing
    forcing[:st_ppm] = co2model[:st_ppm];
    forcing[:F_CH4] = CH4model[:F_CH4];
    forcing[:F_EX] = RCPmodel[:F_EX];

    # Setup temperature model
    temperaturemodel[:F] = forcing[:F];
    posttemp[:T_AT] = temperaturemodel[:T_AT];
    posttemp[:co2_cum] = co2model[:co2_cum];

    # Setup SLR model
    slr[:T_AT] = temperaturemodel[:T_AT];

    # Setup pattern scaling
    pattscale[:T_AT] = temperaturemodel[:T_AT];

    # Setup Consumption
    cons[:T_country] = pattscale[:T_country];
    cons[:SLR] = slr[:SLR];

    # Setup Utility
    utility[:conspc] = cons[:conspc];

    model
end

function full_model()
    model = test_model();

    RCPmodel = addRCP(model, "RCP8.5");
    co2model = addCO2Model(model, "Expectation");
    CH4model = addCH4Model(model, "Value");
    forcing = addForcing(model, "Best fit multi-model mean");
    safmodel = addSAFModel(model, "Distribution mean");
    temperaturemodel = addTemperatureModel(model, "Best fit multi-model mean");
    posttemp = addPostTemperature(model, "Expectation");
    interact = addInteractions(model);
    pcfmodel = addPCFModel(model, "Fit of Hope and Schaefer (2016)");
    omh = addOMH(model, "Whiteman et al. beta 20 years");
    amaz = addAmazonDieback(model, "Cai et al. central value");
    gismodel = addGISModel(model, "Nordhaus central value");
    wais = addWAISmodel(model, "Value");
    slr = addSLR(model);
    ismmodel = addISMModel(model, "Value");
    pattscale = addPatternScaling(model);
    amoc = addAMOC(model, "Hadley");
    cons = addConsumption(model, "none", "none", "SSP2");
    nonmarket = addNonMarketDamages(model);
    utility = addUtility(model, "SSP2");

    # Setup CO2 model
    co2model[:co2_rcp] = RCPmodel[:co2_rcp];
    co2model[:co2_pcf] = pcfmodel[:CO2_PF];
    co2model[:co2_amazon] = amaz[:CO2_AMAZ];
    co2model[:alpha] = posttemp[:alpha];

    # Setup CH4 model
    CH4model[:ch4_rcp] = RCPmodel[:ch4_rcp];
    CH4model[:ch4_pcf] = pcfmodel[:CH4_PF];
    CH4model[:ch4_omh] = omh[:CH4_OMH];
    CH4model[:ch4_conc_rcp] = RCPmodel[:ch4_conc_rcp];
    CH4model[:n2o_conc_rcp] = RCPmodel[:n2o_conc_rcp];

    # Setup forcing
    forcing[:st_ppm] = co2model[:st_ppm];
    forcing[:F_CH4] = CH4model[:F_CH4];
    forcing[:F_EX] = RCPmodel[:F_EX];

    # Setup SAF model
    safmodel[:F] = forcing[:F];

    # Setup temperature model
    temperaturemodel[:F] = forcing[:F];
    temperaturemodel[:T_AT_adjustment] = safmodel[:T_AT_adjustment];
    posttemp[:T_AT] = temperaturemodel[:T_AT];
    posttemp[:co2_cum] = co2model[:co2_cum];

    # Setup Interactions
    interact[:I_AMOC] = amoc[:I_AMOC];
    interact[:VGIS] = gismodel[:VGIS];
    interact[:p_WAIS] = wais[:p_WAIS];
    interact[:I_AMAZ] = amaz[:I_AMAZ];
    interact[:mNINO3pt4] = ismmodel[:mNINO3pt4];

    # Setup PCF model
    pcfmodel[:T_AT] = temperaturemodel[:T_AT];

    # Setup OMH model
    omh[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
    omh[:T_AT] = temperaturemodel[:T_AT];

    # Setup WAIS model
    wais[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
    wais[:T_AT] = temperaturemodel[:T_AT];
    wais[:f_WAIS] = interact[:f_WAIS];

    # Setup Amazon Dieback model
    amaz[:T_AT] = temperaturemodel[:T_AT];
    amaz[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
    amaz[:probmult] = interact[:f_AMAZ];

    # Setup GIS model
    gismodel[:T_AT] = temperaturemodel[:T_AT];
    gismodel[:f_GIS] = interact[:f_GIS];

    # Setup SLR model
    slr[:T_AT] = temperaturemodel[:T_AT];

    # Setup ISM model
    ismmodel[:T_AT] = temperaturemodel[:T_AT];
    ismmodel[:st_ppm] = co2model[:st_ppm];
    ismmodel[:uniforms] = reshape(rand(Uniform(0, 1), dim_count(model, :time) * dim_count(model, :monsoonsteps)),
                                  dim_count(model, :time), dim_count(model, :monsoonsteps));
    ismmodel[:SO_2] = RCPmodel[:SO_2];
    ismmodel[:f_NINO] = interact[:f_NINO];
    pattscale[:T_AT] = temperaturemodel[:T_AT];

    # Setup AMOC
    amoc[:T_AT] = temperaturemodel[:T_AT];
    amoc[:scale_country] = pattscale[:T_country];
    amoc[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
    amoc[:f_AMOC] = interact[:f_AMOC];

    # Setup Consumption
    cons[:T_country] = pattscale[:T_country];
    cons[:SLR] = slr[:SLR];

    # Setup Non-market damages
    nonmarket[:T_AT] = temperaturemodel[:T_AT];
    nonmarket[:conspc] = cons[:conspc];

    # Setup Utility
    utility[:conspc] = cons[:conspc];
    utility[:lossfactor] = nonmarket[:lossfactor];

    model
end

model = base_model()
run(model)

model = full_model()
run(model)
