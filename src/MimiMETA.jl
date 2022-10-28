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

function base_model(; rcp="RCP4.5", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="none", slrdamage="none")
    model = test_model();

    rcpmodel = addRCP(model, rcp);
    co2model = addCO2Model(model, co2);
    ch4model = addCH4Model(model, ch4);
    forcing = addForcing(model, warming);
    temperaturemodel = addTemperatureModel(model, warming);
    posttemp = addPostTemperature(model, co2);
    slr = addSLR(model);
    pattscale = addPatternScaling(model);
    cons = addConsumption(model, tdamage, slrdamage, ssp);
    utility = addUtility(model, ssp);

    # Setup CO2 model
    co2model[:co2_rcp] = rcpmodel[:co2_rcp];
    co2model[:co2_alpha] = posttemp[:alpha];

    # Setup CH4 model
    ch4model[:ch4_rcp] = rcpmodel[:ch4_rcp];
    ch4model[:ch4_conc_rcp] = rcpmodel[:ch4_conc_rcp];
    ch4model[:n2o_conc_rcp] = rcpmodel[:n2o_conc_rcp];

    # Setup forcing
    forcing[:st_ppm] = co2model[:st_ppm];
    forcing[:F_CH4] = ch4model[:F_CH4];
    forcing[:F_EX] = rcpmodel[:F_EX];

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

function full_model(; rcp="RCP4.5", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="pointestimate", slrdamage="mode", saf="Distribution mean", interaction=true, pcf="Fit of Hope and Schaefer (2016)", omh="Whiteman et al. beta 20 years", amaz="Cai et al. central value", gis="Nordhaus central value", wais="Value", ism="Value", amoc="IPSL", nonmarketdamage=false)
    model = base_model(rcp=rcp, ssp=ssp, co2=co2, ch4=ch4, warming=warming, tdamage=tdamage, slrdamage=slrdamage);

    if saf != false
        safmodel = addSAFModel(model, saf, before=:TemperatureModel);

        connect_param!(model, :SAFModel=>:F, :Forcing=>:F);
        connect_param!(model, :TemperatureModel=>:T_AT_adjustment, :SAFModel=>:T_AT_adjustment);
    end
    if interaction != false
        interact = addInteractions(model, after=:PostTemperature);
    end
    if pcf != false
        pcfmodel = addPCFModel(model, pcf, after=:TemperatureModel);

        connect_param!(model, :CO2Model=>:co2_pcf, :PCFModel=>:CO2_PF);
        connect_param!(model, :CH4Model=>:ch4_pcf, :PCFModel=>:CH4_PF);

        connect_param!(model, :PCFModel=>:T_AT, :TemperatureModel=>:T_AT);
    end
    if omh != false
        omhmodel = addOMH(model, omh, after=:TemperatureModel);

        connect_param!(model, :CH4Model=>:ch4_omh, :OMH=>:CH4_OMH);

        omhmodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
        connect_param!(model, :OMH=>:T_AT, :TemperatureModel=>:T_AT);
    end
    if amaz != false
        amazmodel = addAmazonDieback(model, amaz, after=ifelse(interaction, :Interactions, :TemperatureModel));

        connect_param!(model, :CO2Model=>:co2_amazon, :AmazonDieback=>:CO2_AMAZ);

        connect_param!(model, :AmazonDieback=>:T_AT, :TemperatureModel=>:T_AT);
        amazmodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
        if interaction != false
            amazmodel[:probmult] = interact[:f_AMAZ];
        end
    end
    if gis != false
        gismodel = addGISModel(model, gis, after=ifelse(interaction, :Interactions, :TemperatureModel));

        connect_param!(model, :GISModel=>:T_AT, :TemperatureModel=>:T_AT);
        if interaction != false
            gismodel[:f_GIS] = interact[:f_GIS];
        end
    end
    if wais != false
        waismodel = addWAISmodel(model, wais, after=ifelse(interaction, :Interactions, :TemperatureModel));

        waismodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
        connect_param!(model, :WAISmodel=>:T_AT, :TemperatureModel=>:T_AT);
        if interaction != false
            waismodel[:f_WAIS] = interact[:f_WAIS];
        end
    end
    if ism != false
        ismmodel = addISMModel(model, ism, after=ifelse(interaction, :Interactions, :TemperatureModel));

        connect_param!(model, :ISMModel=>:T_AT, :TemperatureModel=>:T_AT);
        connect_param!(model, :ISMModel=>:st_ppm, :CO2Model=>:st_ppm);
        ismmodel[:uniforms] = reshape(rand(Uniform(0, 1), dim_count(model, :time) * dim_count(model, :monsoonsteps)),
                                      dim_count(model, :time), dim_count(model, :monsoonsteps));
        connect_param!(model, :ISMModel=>:SO_2, :RCP=>:SO_2);
        if interaction != false
            ismmodel[:f_NINO] = interact[:f_NINO];
        end
    end
    if amoc != false
        amocmodel = addAMOC(model, amoc, after=:PatternScaling);

        connect_param!(model, :AMOC=>:T_AT, :TemperatureModel=>:T_AT);
        connect_param!(model, :AMOC=>:scale_country, :PatternScaling=>:T_country);
        amocmodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
        if interaction != false
            amocmodel[:f_AMOC] = interact[:f_AMOC];
        end
    end
    if nonmarketdamage != false
        nonmarket = addNonMarketDamages(model, after=:Consumption);

        connect_param!(model, :NonMarketDamages=>:T_AT, :TemperatureModel=>:T_AT);
        connect_param!(model, :NonMarketDamages=>:conspc, :Consumption=>:conspc);

        connect_param!(model, :Utility=>:lossfactor, :NonMarketDamages=>:lossfactor);
    end

    if interaction != false
        # Setup Interactions
        if amoc != false
            interact[:I_AMOC] = amocmodel[:I_AMOC];
        end
        if gis != false
            interact[:VGIS] = gismodel[:VGIS];
        end
        if wais != false
            interact[:p_WAIS] = waismodel[:p_WAIS];
        end
        if amaz != false
            interact[:I_AMAZ] = amazmodel[:I_AMAZ];
        end
        if ism != false
            interact[:mNINO3pt4] = ismmodel[:mNINO3pt4];
        end
    end

    model
end

## TODO: Need to move these into a unit test
# model = base_model()
# run(model)

# model = full_model()
# run(model)
