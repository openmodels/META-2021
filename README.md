# META 2021
The Model for Economic Tipping point Analysis

META 2021 is an advanced integrated assessment model (SC-IAM), designed as a model-based meta-analysis of the effects of tipping points on the social cost of carbon (SCC). The model simulates greenhouse gas emissions, temperature and sea-level rise, and market and non-market damages at the country-level, and the effects of eight climate tipping points that have been studied in the climate economics literature.

META 2021 is introduced in; **Dietz, Rising, Stoerk, and Wagner (2021): "Economic impacts of tipping points in the climate system", PNAS.** [LINK TO DOI]

See that paper and its supplementary information for further details. Please cite the paper when using META in your research.

## Installation and Use

The main model is implemented as an Excel spreadsheet, which uses [@RISK](https://www.palisade.com/risk/) for Monte Carlo analysis. You should have @RISK installed to use META 2021.

To open the model, first load @RISK, and within @RISK load `Model/META model July 2021.xlsx`.

The model links to a number of supporting files. You should therefore also open the following files:

 - `BHM betas v3pd.xlsx`: Country-level damage function coefficients, derived from [Burke et al. (2015)](https://www.nature.com/articles/nature15725).
 - `National gdp per capita ppp.xlsx`: Projected GDP per capita from the Shared Socioeconomic Pathways (SSP).
 - `National population.xlsx`: Projected population totals from the Shared Socioeconomic Pathways (SSP).

The first tab, `Settings`, allows you to activate tipping points, choose emissions and socioeconomic scenarios, and choose other settings. It also contains a dashboard, which provides an overview of your settings, including key parameters and sets of parameters (see below).

Input parameters and their uncertainty are defined in the `Parameters` tab. Cells that are intended to be changed, in order to switch from one parameter scheme to another, are shaded yellow.

All other tabs compute intermediate and final results, including global temperature rise (`Temperatures` column G), sea-level rise (`SLR` column D), country-level downscaled temperature (`Pattern scaling`), income levels after impacts (`National cons per cap`), and global social welfare (`Welfare & SCCO2 calculator` cell C2).

## Contents of the repository

All model files are in the `Model` folder, and consist of:

 - `META model July 2021.xlsx`: The main model file (see Installation and Use).
 - `API_NY.GDP.PCAP.PP.CD_DS2_en_excel_v2_126211.xls`: World Bank historical GDP per capita.
 - `API_NY.GNS.ICTR.ZS_DS2_en_excel_v2_41442.xls`: World Bank historical savings rates.
 - `API_SP.POP.TOTL_DS2_en_excel_v2_126122.xls`: World Bank historical population totals.
 - `BHM betas v3pd.xlsx`: Country-level damage function coefficients, derived from [Burke et al. (2015)](https://www.nature.com/articles/nature15725).
 - `National gdp per capita ppp.xlsx`: Projected GDP per capita from the Shared Socioeconomic Pathways (SSP).
 - `National population.xlsx`: Projected population totals from the Shared Socioeconomic Pathways (SSP).

The main model includes links to the other files in the `Model` folder.

## Computing the social cost of carbon

To calculate the SCC, run the model without the additional pulse of CO2 (set `Settings` cell C14 to `No`) and record welfare (`Welfare & SCCO2 calculator` cell C2). Then run the model again with the additional pulse of CO2 (set `Settings` cell C14 to `Yes`) and record welfare (`Welfare & SCCO2 calculator` cell C2 again). Calculate the difference in welfare. Then divide by ∂W/∂C(2020), using world mean consumption/capita in `Welfare & SCCO2 calculator` cell C3. You may wish to inflate the resulting number to current dollars!

## Running the model without @RISK

We very strongly recommend installing @RISK before loading META 2021. However, there is a setting that enables the models' many random parameters to be treated as deterministic (set `Settings` cell C16 to `No`). The model should then work in standard Excel but will not perform Monte Carlo simulations.

## Model details

### Settings

![image](https://user-images.githubusercontent.com/579448/125863417-3a22da1a-2371-46c2-81ad-04a12c1182ae.png)

The top of the `Settings` tab provides a dashboard of high-level choices:

 - Each tipping point (under `B2`) can be independently turned `Off` or `On`. Each tipping point is discussed in more detail below. The `Interactions` option specifies if the triggering of one tipping point changes the probability of other tipping points.
 - The `SCC` section (under `B13`) has options for including an extra pulse of emissions (`Extra tonne`), including uncertainty (`Stochastic`), and calculating damages from climate change (`Climate impacts`).
 - The `RCP scenario` (at `F3`) sets the emissions scenario from amongst the CMIP5 representative concentration pathways.
 - The `SSP scenario` (at `F9`) sets the socioeconomic scenario from amongst the Shared Socioeconomic Pathways.
 - `Non-market damages` can be calculated, in addition to the standard Burke et al. market damages, if `F3` is set to `MERGE`.
 - Four different representations of the AMOC tipping point are included, as selected in cell `L3`.

The bottom of the sheet draws from the `Parameters` sheet, based on the selections above.

### Parameters

![image](https://user-images.githubusercontent.com/579448/125862653-d9eae9e6-8977-4903-94d4-d80cec2a645b.png)

Column B contains inputs to the model, columns C onwards to the right provide various alternatives that you can choose. So, you set the yellow cells in column B equal to whichever input alternative (i.e. column) you want to use. Where available, 'distribution' gives you the stochastic option. Also change the title cell (e.g. cell B4 for the carbon cycle module) so that the dashboard is updated.
