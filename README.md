# META 2021
The Model for Economic Tipping (point) Analysis

META 2021 is an advanced integrated assessment model (SC-IAM),
designed as a model-based meta-analysis of the effects of tipping
points on the social cost of carbon (SCC). The model simulates
greenhouse gas emissions, temperature and sea-level rise, market and
non-market damages at the country-level, and the effects of eight
climate tipping points that have been implemented in the climate
economics literature.

META 2021 is introduced in **Dietz, Rising, Stoerk, and Wagner (2021): "Economic impacts of tipping points in the climate system", PNAS.** [LINK TO DOI]

See that paper and its supplemental information for more information. Please cite that paper when using META in your research.

## Installation and Use

The main model is implemented as an Excel spreadsheet, which uses
[@RISK](https://www.palisade.com/risk/) for Monte Carlo analysis. You
must have @RISK installed to use META 2021.

To open the model, first load @RISK, and within @RISK load `Model/META
model July 2021.xlsx`.

The first tab, `Settings`, provides a dashboard for selecting tipping
points, emissions and socioeconomic scenarios, and other settings.

Input parameters and their uncertainty are defined on the `Parameters`
tab.

All other tabs compute intermediate and final results, including
global temperature rise (`Temperatures` column G), sea-level rise
(`SLR` column D), country-level downscaled temperature (`Pattern
scaling`), income levels after impacts (`National cons per cap`), and
global social welfare (`Welfare & SCCO2 calculator` cell C2).

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
