using Distributions, CSVFiles
using MimiFAIRv2

include("../src/basemodel.jl")
include("../src/components/RCP.jl")
include("../src/components/CO2Converter.jl")
include("../src/components/CH4Converter.jl")
include("../src/components/TemperatureConverter.jl")
include("../src/components/SAF.jl")
include("../src/components/temperature_withadj.jl")
include("../src/components/PCF.jl")
include("../src/components/OMH.jl")
include("../src/components/AmazonDieback.jl")
include("../src/components/GIS.jl")
include("../src/components/WAIS.jl")
include("../src/components/AIS.jl")
include("../src/components/SLR.jl")
include("../src/components/ISM.jl")
include("../src/components/PatternScaling.jl")
include("../src/components/AMOC.jl")
include("../src/components/Interactions.jl")
include("../src/components/Consumption.jl")
include("../src/components/NonMarketDamages.jl")
include("../src/components/Utility.jl")
include("../src/components/TotalDamages.jl")
include("../src/components/BGE.jl")
#include("../src/components/DEBUG.jl")

function base_model(; rcp="CP-Base", ssp="SSP2", tdamage="none", slrdamage="none")
    model = MimiFAIRv2.get_model(end_year=2200)

    # Get backup data
    emissions_data = DataFrame(load(joinpath(dirname(pathof(MimiFAIRv2)), "..", "data", "rcmip_ssp245_emissions_1750_to_2500.csv"), skiplines_begin=6))

    set_special_model_dimensions!(model)

    rcpmodel = addRCP(model, rcp, before=:co2_cycle);
    co2converter = addCO2Converter(model, after=:RCP);
    ch4converter = addCH4Converter(model, after=:RCP);
    tempconverter = addTemperatureConverter(model, after=:temperature);
    slr = addSLR(model);
    pattscale = addPatternScaling(model);
    cons = addConsumption(model, tdamage, slrdamage, ssp);
    utility = addUtility(model, ssp);
    damages = addTotalDamages(model);
    bge_comp = addBGE(model)
    #debug = addDEBUG(model)

    # Setup CO2 converter
    co2converter[:co2_rcp] = rcpmodel[:co2_rcp];

    # Setup CH4 converter
    ch4converter[:ch4_rcp] = rcpmodel[:ch4_rcp];
    ch4converter[:ch4_2009] = emissions_data.methane[emissions_data[!, 1] .== 2009][1]

    # Setup temperature converter
    connect_param!(model, :TemperatureConverter => :T, :temperature => :T)

    # Feed converters into FAIR
    connect_param!(model, :co2_cycle => :E_co2, :CO2Converter => :E_co2, emissions_data.carbon_dioxide[emissions_data[!, 1] .<= 2200])
    connect_param!(model, :ch4_cycle => :E_ch4, :CH4Converter => :E_ch4, emissions_data.methane[emissions_data[!, 1] .<= 2200])

    # Setup SLR model
    connect_param!(model, :SLRModel => :T_AT, :TemperatureConverter => :T_AT)

    # Setup pattern scaling
    connect_param!(model, :PatternScaling => :T_AT, :TemperatureConverter => :T_AT)

    # Setup Consumption
    cons[:T_country] = pattscale[:T_country];
    cons[:SLR] = slr[:SLR];

    # Setup Utility
    utility[:conspc] = cons[:conspc];

    # Setup TotalDamages
    damages[:population] = utility[:pop];
    damages[:postdamage_consumption_percap_percountry] = cons[:conspc];
    damages[:baseline_consumption_percap_percountry] = cons[:baseline_consumption_percap_percountry];

    # Setup BGE
    bge_comp[:pop] = utility[:pop];
    bge_comp[:utility] = utility[:utility];
    bge_comp[:baseline_consumption_percap_percountry] = cons[:baseline_consumption_percap_percountry];

    model
end

function full_model(; rcp="CP-Base", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="pointestimate", slrdamage="mode", saf="Distribution mean", interaction=true, pcf="Fit of Hope and Schaefer (2016)", omh="Whiteman et al. beta 20 years", amaz="Cai et al. central value", gis="Nordhaus central value", ais="AIS", ism="Value", amoc="IPSL", nonmarketdamage=false)
    model = base_model(rcp=rcp, ssp=ssp, tdamage=tdamage, slrdamage=slrdamage);

    if saf != false
        safmodel = addSAFModel(model, saf, before=:temperature);

        replace!(model, :temperature => temperature_withadj)

        connect_param!(model, :SAFModel=>:F, :radiative_forcing => :total_RF);
        connect_param!(model, :temperature=>:T1_adjustment, :SAFModel=>:T_AT_adjustment, zeros(2200 - 1750 + 1), ignoreunits=true);
    end
    if interaction != false
        interact = addInteractions(model, after=:TemperatureConverter);
    end
    if pcf != false
        pcfmodel = addPCFModel(model, pcf, after=:TemperatureConverter);

        connect_param!(model, :CO2Converter=>:co2_pcf, :PCFModel=>:CO2_PF);
        connect_param!(model, :CH4Converter=>:ch4_pcf, :PCFModel=>:CH4_PF);

        connect_param!(model, :PCFModel=>:T_AT, :TemperatureConverter => :T_AT);
    end
    if omh != false
        omhmodel = addOMH(model, omh, after=:TemperatureConverter);

        connect_param!(model, :CH4Converter=>:ch4_omh, :OMH=>:CH4_OMH);

        omhmodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
        connect_param!(model, :OMH=>:T_AT, :TemperatureConverter => :T_AT);
    end
    if amaz != false
        amazmodel = addAmazonDieback(model, amaz, after=ifelse(interaction, :Interactions, :TemperatureConverter));

        connect_param!(model, :CO2Converter=>:co2_amazon, :AmazonDieback=>:CO2_AMAZ);

        connect_param!(model, :AmazonDieback=>:T_AT, :TemperatureConverter => :T_AT);
        amazmodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
        if interaction != false
            amazmodel[:probmult] = interact[:f_AMAZ];
        end
    end
    if gis != false
        gismodel = addGISModel(model, gis, after=ifelse(interaction, :Interactions, :TemperatureConverter));

        connect_param!(model, :GISModel=>:T_AT, :TemperatureConverter => :T_AT);
        if interaction != false
            gismodel[:f_GIS] = interact[:f_GIS];
        end
        connect_param!(model, :SLRModel=>:SLR_GIS, :GISModel=>:SLR_GIS);
    end
    if ais == "AIS"
        aismodel = addAISmodel(model, after=ifelse(interaction, :Interactions, :TemperatureConverter))
        connect_param!(model, :AISmodel=>:T_AT, :TemperatureConverter => :T_AT);
        connect_param!(model, :AISmodel=>:T_AT_tminus100, :TemperatureConverter => :T_AT_tminus100);
        if interaction != false
            aismodel[:f_AIS] = interact[:f_AIS];
        end
        connect_param!(model, :SLRModel=>:SLR_AIS, :AISmodel=>:SLR_AIS);
    elseif ais == "WAIS"
        waismodel = addWAISmodel(model, after=ifelse(interaction, :Interactions, :TemperatureConverter));
        waismodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
        connect_param!(model, :WAISmodel=>:T_AT, :TemperatureConverter => :T_AT);
        if interaction != false
            waismodel[:f_WAIS] = interact[:f_WAIS];
        end
        connect_param!(model, :SLRModel=>:SLR_AIS, :WAISmodel=>:SLR_WAIS);
    end
    if ism != false
        ismmodel = addISMModel(model, ism, after=ifelse(interaction, :Interactions, :TemperatureConverter));

        connect_param!(model, :ISMModel=>:T_AT, :TemperatureConverter => :T_AT);
        connect_param!(model, :ISMModel=>:st_ppm, :co2_cycle=>:co2);
        ismmodel[:uniforms] = reshape(rand(Uniform(0, 1), dim_count(model, :time) * dim_count(model, :monsoonsteps)),
                                      dim_count(model, :time), dim_count(model, :monsoonsteps));
        connect_param!(model, :ISMModel=>:SO_2, :RCP=>:SO_2);
        if interaction != false
            ismmodel[:f_NINO] = interact[:f_NINO];
        end
        connect_param!(model, :Consumption=>:extradamage, :ISMModel=>:extradamage);
    end
    if amoc != false
        amocmodel = addAMOC(model, amoc, after=:PatternScaling);

        connect_param!(model, :AMOC=>:T_AT, :TemperatureConverter => :T_AT);
        connect_param!(model, :AMOC=>:scale_country, :PatternScaling=>:T_country);
        amocmodel[:uniforms] = rand(Uniform(0, 1), dim_count(model, :time));
        if interaction != false
            amocmodel[:f_AMOC] = interact[:f_AMOC];
        end
    end
    if nonmarketdamage != false
        nonmarket = addNonMarketDamages(model, after=:Consumption);

        connect_param!(model, :NonMarketDamages=>:T_AT, :TemperatureConverter => :T_AT);
        connect_param!(model, :NonMarketDamages=>:conspc, :Consumption=>:conspc);

        connect_param!(model, :Utility=>:lossfactor, :NonMarketDamages=>:lossfactor);
        connect_param!(model, :BGE=>:lossfactor, :NonMarketDamages=>:lossfactor);
    end

    if interaction != false
        # Setup Interactions
        if amoc != false
            interact[:I_AMOC] = amocmodel[:I_AMOC];
        end
        if gis != false
            interact[:VGIS] = gismodel[:VGIS];
        end
        if ais == "AIS"
            interact[:totalSLR_Ross] = aismodel[:totalSLR_Ross];
            interact[:totalSLR_Amundsen] = aismodel[:totalSLR_Amundsen];
            interact[:totalSLR_Weddell] = aismodel[:totalSLR_Weddell];
            interact[:totalSLR_Peninsula] = aismodel[:totalSLR_Peninsula];
        elseif ais == "WAIS"
            interact[:I_WAIS] = waismodel[:I_WAIS];
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
