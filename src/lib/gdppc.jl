using Statistics

## Historical socioeconomics

gdppcs = CSV.read("../data/GDPpc.csv", DataFrame)
pops = CSV.read("../data/Population.csv", DataFrame)

function getgdppc(iso, year)
    index = findfirst(gdppcs.ISO .== iso)

    if isnothing(index)
        missing
    else
        gdppcs[index, Symbol(year)]
    end
end

function getpop(iso, year)
    index = findfirst(pops.ISO .== iso)

    if isnothing(index)
        missing
    else
        pops[index, Symbol(year)]
    end
end

## Projected socioeconomics

gdps_ssp = CSV.read("../data/GDP-SSP.csv", DataFrame)
pops_ssp = CSV.read("../data/Pop-SSP.csv", DataFrame)

@assert all(gdps_ssp.SCENARIO .== pops_ssp.SCENARIO)
@assert all(gdps_ssp.REGION .== pops_ssp.REGION)
@assert all(names(gdps_ssp) .== names(pops_ssp))

gdps_ssp.SCENARIO_SHORT = [scenario[1:4] for scenario in gdps_ssp.SCENARIO]
gdps_ssp.REGION_SHORT = [scenario[5:end] for scenario in gdps_ssp.REGION]

countryinfo = CSV.read("../data/countryinfo.csv", DataFrame)

function getregion(iso)
    index = findfirst(countryinfo.ISO .== iso)
    if isnothing(index)
        missing
    else
        countryinfo.REG[index]
    end
end

function getgdppc_ssp(region, scenario, year)
    if ismissing(region)
        return 0
    end

    if year == 2005 || year in 2010:10:2100
        indices = findall((gdps_ssp.REGION_SHORT .== region) .& (gdps_ssp.SCENARIO_SHORT .== scenario))
        yrsym = Symbol(year)

        mean(1000 * gdps_ssp[indices, yrsym] ./ pops_ssp[indices, yrsym])
    elseif year > 2005 && year < 2010
        baseline = getgdppc_ssp(region, scenario, 2005)
        baseline + (year - 2005) * (getgdppc_ssp(region, scenario, 2010) - baseline) / 5
    elseif year > 2100
        getgdppc_ssp(region, scenario, 2100)
    else
        baseyear = 10 * trunc(Int, year / 10)
        baseline = getgdppc_ssp(region, scenario, baseyear)
        baseline + (year - baseyear) * (getgdppc_ssp(region, scenario, (baseyear + 10)) - baseline) / 10
    end
end

function getpopweight_ssp(region, scenario)
    getpop_ssp(region, scenario, 2015) / sum([getpop_ssp(reg, scenario, 2015) for reg in unique(gdps_ssp.REGION_SHORT)])
end

function getpop_ssp(region, scenario, year)
    if year == 2005 || year in 2010:10:2100
        yrsym = Symbol(year)
        indices = findall((gdps_ssp.REGION_SHORT .== region) .& (gdps_ssp.SCENARIO_SHORT .== scenario) .& (pops_ssp[!, yrsym] .!= 0.0))

        mean(pops_ssp[indices, yrsym])
    elseif year < 2005
        error("Population years before 2005 or after 2100 are not implemented.")
    elseif year > 2100
        getpop_ssp(region, scenario, 2100)
    else
        baseyear = 10 * trunc(Int, year / 10)
        baseline = getpop_ssp(region, scenario, baseyear)
        baseline + (year - baseyear) * (getpop_ssp(region, scenario, (baseyear + 10)) - baseline) / 10
    end
end
