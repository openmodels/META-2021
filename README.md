# META 2021
The Model for Economic Tipping point Analysis

META 2021 is an advanced integrated assessment model (SC-IAM), designed as a model-based meta-analysis of the effects of tipping points on the social cost of carbon (SCC). The model simulates greenhouse gas emissions, temperature and sea-level rise, and market and non-market damages at the country level, and the effects of eight climate tipping points that have been studied in the climate economics literature.

META 2021 is introduced in; **Dietz, Rising, Stoerk, and Wagner (2021): "Economic impacts of tipping points in the climate system", PNAS, 118(34), e2103081118.** [https://doi.org/10.1073/pnas.2103081118]

See that paper and its supplementary information for further
details. Please cite the paper when using META in your research.

## Model versions

The repository contains two versions of the model, which produce
identical results to high precision. The initial version of the model
is implemented in Excel using the @RISK extension to perform Monte
Carlo simulations. We have now produced a version implemented in Mimi
[https://www.mimiframework.org/], an integrated assessment modeling
framework developed in Julia [https://julialang.org/].

The Excel version is more user-friendly, but additional analyses may
require the Mimi version. The Mimi version will be used for future
developments.

## Excel @RISK model

### Installation and Use

It is recommended that you download this entire repository to use META, since the model relies on a number of linked files. You can do this at the `Code` link above or download a [zip archive](https://github.com/openmodels/META-2021/archive/refs/heads/master.zip).

The main model is implemented as an Excel spreadsheet, which uses [@RISK](https://www.palisade.com/risk/) for Monte Carlo analysis. You should have @RISK installed to use META 2021.

To open the model, first load @RISK, and within @RISK load `excel/META model July 2021.xlsx`.

The model links to a number of supporting files. You should therefore also open the following files:

 - `BHM betas v3pd.xlsx`: Country-level damage function coefficients, derived from [Burke et al. (2015)](https://www.nature.com/articles/nature15725).
 - `National gdp per capita ppp.xlsx`: Projected GDP per capita from the Shared Socioeconomic Pathways (SSPs).
 - `National population.xlsx`: Projected population totals from the Shared Socioeconomic Pathways (SSPs).

The first tab, `Settings`, allows you to activate tipping points, choose emissions and socioeconomic scenarios, and choose other settings. It also contains a dashboard, which provides an overview of your settings, including key parameters and sets of parameters (see below).

Input parameters and their uncertainty are defined in the `Parameters` tab. Cells that are intended to be changed, in order to switch from one parameter scheme to another, are shaded yellow.

All other tabs compute intermediate and final results, including global temperature rise (`Temperatures` column G), sea-level rise (`SLR` column D), country-level downscaled temperature (`Pattern scaling`), income levels after impacts (`National cons per cap`), and global social welfare (`Welfare & SCCO2 calculator` cell C2).

### Contents of the Excel model

All model files are in the `excel` folder, and consist of:

 - `META model July 2021.xlsx`: The main model file (see Installation and Use).
 - `API_NY.GDP.PCAP.PP.CD_DS2_en_excel_v2_126211.xls`: World Bank historical GDP per capita.
 - `API_NY.GNS.ICTR.ZS_DS2_en_excel_v2_41442.xls`: World Bank historical savings rates.
 - `API_SP.POP.TOTL_DS2_en_excel_v2_126122.xls`: World Bank historical population totals.
 - `BHM betas v3pd.xlsx`: Country-level damage function coefficients, derived from [Burke et al. (2015)](https://www.nature.com/articles/nature15725).
 - `National gdp per capita ppp.xlsx`: Projected GDP per capita from the Shared Socioeconomic Pathways (SSPs).
 - `National population.xlsx`: Projected population totals from the Shared Socioeconomic Pathways (SSPs).

The main model includes links to the other files in the `excel` folder.

### Computing the social cost of carbon

To calculate the SCC, run the model without the additional pulse of CO2 (set `Settings` cell C14 to `No`) and record welfare (`Welfare & SCCO2 calculator` cell C2). Then run the model again with the additional pulse of CO2 (set `Settings` cell C14 to `Yes`) and record welfare (`Welfare & SCCO2 calculator` cell C2 again). Calculate the difference in welfare. Then divide by ∂W/∂C(2020), using world mean consumption/capita in `Welfare & SCCO2 calculator` cell C3. You may wish to inflate the resulting number to current dollars!

### Running the model without @RISK

We very strongly recommend installing @RISK before loading META 2021. However, there is a setting that enables the models' many random parameters to be treated as deterministic (set `Settings` cell C16 to `No`). The model should then work in standard Excel but will not perform Monte Carlo simulations.

### Model details

#### Settings

![image](https://user-images.githubusercontent.com/579448/125863417-3a22da1a-2371-46c2-81ad-04a12c1182ae.png)

The top of the `Settings` tab provides a dashboard of high-level choices:

 - Each tipping point (under `B2`) can be independently turned `Off` or `On`. Each tipping point is discussed in more detail below. The `Interactions` option specifies if the triggering of one tipping point changes the probability of other tipping points.
 - The `SCC` section (under `B13`) has options for including an extra pulse of emissions (`Extra tonne`), including uncertainty (`Stochastic`), and calculating damages from climate change (`Climate impacts`).
 - The `RCP scenario` (at `F3`) sets the emissions scenario from amongst the CMIP5 representative concentration pathways.
 - The `SSP scenario` (at `F9`) sets the socioeconomic scenario from amongst the Shared Socioeconomic Pathways.
 - `Non-market damages` can be calculated, in addition to the standard Burke et al. market damages, if `F3` is set to `MERGE`.
 - Four different representations of the AMOC tipping point are included, as selected in cell `L3`.

The bottom of the sheet draws from the `Parameters` sheet, based on the selections above.

#### Parameters

![image](https://user-images.githubusercontent.com/579448/125862653-d9eae9e6-8977-4903-94d4-d80cec2a645b.png)

Column B contains inputs to the model, columns C onwards to the right provide various alternatives that you can choose. So, you set the yellow cells in column B equal to whichever input alternative (i.e. column) you want to use. Where available, 'distribution' gives you the stochastic option. Also change the title cell (e.g. cell B4 for the carbon cycle module) so that the dashboard is updated.

## Mimi model

See the description of the Excel model for details on the
component-based structure and parameterization, which are maintained
in the Mimi model.

### Directories in the repository

The following directories are used for the Mimi model:
 - `data`: Input parameters and validation outputs (under
   `data/benchmark`).
 - `src`: The model code. The scripts directly contained in this directory
   support various types of analyses, with internal model code in
   subdirectories. Specifically, the model components are under
   `src/components` and additional functions are in `src/lib`.
 - `test`: Unit tests, for each component, for the system-wide
   results, and for the Monte Carlo system.
   
 ### Basic use cases
 
 Please note that all code is designed to be run with the working
 directory set to a subdirectory of the repository (e.g., `src` or you
 can create a subdirectory `analysis`).
 
 #### 1. Running the full deterministic model
 
 The full model is constructed using `full_model(...)`, defined in
 `src/MimiMETA.jl`. The `full_model` function can be called with no
 arguments, to use the default construction, or override the defaults
 with the following arguments:
 
  - `rcp`: Emissions scenario; one of RCP3-PD/2.6, RCP4.5 (default), RCP6, or RCP8.5.
  - `ssp`: Socioeconomic scenario; one of SSP1, SSP2 (default), SSP3, SSP4, SSP5.
  - `co2`: CO2 model calibration; one of AR5-IR, AR5-PI, MESMO (lowest
    decay), ACC2 (highest decay), Expectation (default), or
    Distribution.
  - `ch4`: CH4 model calibration; one of default (default), low, or
    high.
  - `warming`: Forcing model calibration: one of Best fit multi-model
    mean (default), HadGEM2-ES (hottest model), or GISS-E2-R (coldest
    model).
  - `tdamage`: Temperature damages; one of none, distribution,
    pointestimate (default), low, or high.
  - `slrdamage`: Sea-level rise damages; one of none, distribution,
    mode (default), low, or high.
  - `nonmarketdamage`: Non-market damages; May be false (to not use,
    default) or true.
  - `saf`: Surface albedo feedback calibration; May be false (to not
    use) or Distribution mean (default)
  - `pcf`: Permafrost carbon feedback calibration; May be false (to not
    use) or one of Fit of Hope and Schaefer (2016), Kessler central
    value, Kessler 2.5%, Kessler 97.5%, Fit of Hope and Schaefer (2016)
    (default), or Fit of Yumashev et al. (2019).
  - `omh`: Ocean methane hydrates calibration; May be false (to not
    use) or one of Whiteman et al. beta 20 years (default), Whiteman
    et al. uniform 10 years, "Whiteman et al. triangular, mode 10%, 10
    years", Whiteman et al. beta 10 years, "Ceronsky et al. (2011),
    1.784GtCH4 per year, beta", "Ceronsky et al. (2011), 7.8GtCH4 per
    year, beta", "Ceronsky et al. (2011), 0.2GtCH4 per year, beta",
    Whiteman et al. beta 20 years, Whiteman et al. beta 30 years,
    Whiteman et al. uniform 20 years, or "Whiteman et al. triangular,
    mode 10%, 20 years".
  - `amaz`: Amazon dieback calibration; May be false (to not use) or
    one of Cai et al. central value (default), Cai et al. long, or Cai
    et al. short.
  - `gis`: Greenland icesheet calibration; May be false (to not use)
    or one of Nordhaus central value (default), Robinson, Non-linear
    equilibrium function, Ice/SLR low, or Ice/SLR high.
  - `wais`: West Antarctic icesheet calibration; May be false (to not
    use) or true (default).
  - `ism`: Indian summer monsoon calibration; May be false (to not
    use) or Value (default).
  - `amoc`: Atlantic meridional overturning circulation; May be false
    (to not use) or one of Hadley, BCM, IPSL (default), or HADCM.
  - `interaction`: Tipping point interactions; May be false (to not
    use) or true (default).

There is also a `base_model` function which includes only the
non-tipping-point calibration options.

A basic usage is as follows:

```
include("../src/MimiMETA.jl")
model = full_model(rcp="RCP4.5", ssp="SSP2")
run(model)
explore(model)
```

Other examples are shown in `test/test_system_tp.jl` (for the full
model) and `test/test_system_notp.jl` (for the no-tipping-point
model).

#### 2. Running the full Monte Carlo model

To run the model in Monte Carlo mode, you need to first generate the
simulation parameter values, using the `getsim(...)` function, and
then run the Monte Carlos, using the `runsim(...)` function, both of
which are defined in `src/montecarlo.jl`.

`getsim` takes the following arguments (all are required, and given in
order):
 - `trials`: The number of Monte Carlo simulations.
 - `pcf_calib`: May be "Kessler probabilistic" to draw stochastic
   parameters for the PCF model, or one of the options described in
   the deterministic use case.
 - `amazon_calib`: May be "Distribution" to draw stochastic parameters
   for the Amazon dieback model, or one of the options described in
   the deterministic use case.
 - `gis_calib`: May be "Distribution" to draw stochastic parameters
   for the GIS model, or one of the options described in
   the deterministic use case.
 - `wais_calib`: May be "Distribution" to draw stochastic parameters
   for the WAIS model, or one of the options described in the
   deterministic use case.
 - `saf_calib`: May be "Distribution" to draw stochastic parameters
   for the SAF model, or one of the options described in the
   deterministic use case.
 - `persist_dist`: May be true to draw the level of temperature
   damages persistance stochastically, or false.
 - `emuc_dist`: May be true to draw the level of elasticity of
   marginal utility stochastically, or false.
 - `prtp_dist`: May be true to draw the level of pure rate of time
   preference stochastically, or false.
   
 The `runsim` function takes the following parameters, all of which
 must be provided:
  - `model`: A full Mimi model.
  - `draws`: The result of the `getsim` function.
  - `ism_used`: Set to true if the ISM component is included;
    otherwise false.
  - `omh_used`: Set to true if the OMH component is included;
    otherwise false.
  - `amoc_used`: Set to true if the AMOC component is included;
    otherwise false.
  - `saf_used`: Set to true if the SAF component is included;
    otherwise false.    
  - `amazon_calib`: May be one of the options described in the
   deterministic use case or "none" if the Amazon dieback component is
   excluded.
  - `wais_calib`: Set to "Distribution" if the WAIS component is
   included, or one of the options described in the deterministic use
   case.
  - `save_rvs`: Set to true to save all random variables in the final
    result; otherwise false.

There are also `getmodel_base` and `runsim_base` functions, which
include just the non-tipping-point parameters.

A basic usage is as follows:

```
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")
model = full_model(rcp="RCP4.5", ssp="SSP2")
draws = getsim(500, "Fit of Hope and Schaefer (2016)", # PCF
               "Cai et al. central value", # AMAZ
               "Nordhaus central value", # GIS
               "Distribution", # WAIS
               "Distribution", # SAF
               false, # persit
               false, # emuc
               false) # prtp
results = runsim(model, draws, true, # ism_used
                 true, # omh_used
                 true, # amoc_used
                 true, # saf_used
                 "Cai et al. central value", # AMAZ
                 "Distribution") # WAIS
```

Other examples are included in `test/test_montecarlo_tp.jl`.

#### 3. Calculating the social cost of carbon

The `src/scc.jl` script includes functions that help with the
calculation of the SCC, using the infrastructure within Mimi.

The current standard method is `calculate_scc_mc` which takes the
following parameters (all given, in order):
 - `model`: A version of the META model.
 - `preset_fill`: A function to fill in parameters from a pre-computed
   Monte Carlo collection.
 - `maxrr`: The number of Monte Carlos to perform.
 - `pulse_year`: The year to add an additional pulse of CO2.
 - `pulse_size`: The number of Gt to add.
 - `emuc`: The elasticity of marginal utility to use.
 
A basic usage is:
 
```
include("../src/MimiMETA.jl")
include("../src/lib/presets.jl")
include("../src/scc.jl")
benchmark = CSV.read("../data/benchmark/ExcelMETA-alltp.csv", DataFrame)
model = full_model()
preset_fill(rr) = preset_fill_tp(model, benchmark, rr)
calculate_scc_mc(model, preset_fill, nrow(benchmark), 2020, 10., 1.5) # Runs 500 MC reps.
```
