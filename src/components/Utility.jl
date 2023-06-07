include("../lib/gdppc.jl")

@defcomp Utility begin

    # Variables

    pop_region = Variable(index=[time, region], unit="inhabitants")
    pop_ratio_region = Variable(index=[time, region])
    pop_growth_region = Variable(index=[time, region])
    pop = Variable(index=[time, country], unit="inhabitants") # national population
    world_population = Variable(index=[time], unit="inhabitants")

    utility = Variable(index=[time, country]) # national utility
    utility_sum = Variable(index=[time]) # sum of utility across countries
    world_disc_utility = Variable(index=[time]) # world discounted utility
    equiv_conspc = Variable(index=[time]) # equivalent world consumption per capita

    # Caches
    country2region = Variable{Int64}(index=[country]) # gives region index

    # Parameters

    ssp = Parameter{String}()
    convergerate_pop = Parameter()
    decayrate_pop = Parameter()
    popweights_region = Parameter(index=[region])

    EMUC = Parameter() # elasticity of marginal utility of consumption
    PRTP = Parameter(unit="pct") # pure rate of time preference

    conspc = Parameter(index=[time, country], unit="2010 USD PPP") # national consumption per capita
    lossfactor = Parameter(index=[time, country]) # non-market loss factor

    function run_timestep(pp, vv, dd, tt)
        # Determine populations
        isos = dim_keys(model, :country)

        if is_first(tt)
            for rr in dd.region
                vv.pop_region[tt, rr] = getpop_ssp(dim_keys(model, :region)[rr], pp.ssp, 2010)
                vv.pop_ratio_region[tt, rr] = 1
                vv.pop_growth_region[tt, rr] = 0 # ignored
            end

            for cc in dd.country
                pop = getpop(isos[cc], 2010)
                if ismissing(pop)
                    vv.pop[tt, cc] = 0
                else
                    vv.pop[tt, cc] = pop
                end

                region = getregion(isos[cc])
                vv.country2region[cc] = (ismissing(region) ? 0 : findfirst(dim_keys(model, :region) .== region))
            end
        else
            for rr in dd.region
                vv.pop_region[tt, rr] = getpop_ssp(dim_keys(model, :region)[rr], pp.ssp, gettime(tt))
                vv.pop_ratio_region[tt, rr] = vv.pop_region[tt, rr] / vv.pop_region[TimestepIndex(1), rr]
                if gettime(tt) <= 2100
                    vv.pop_growth_region[tt, rr] = vv.pop_ratio_region[tt, rr] / vv.pop_ratio_region[tt-1, rr] - 1
                else
                    vv.pop_growth_region[tt, rr] = (1-pp.convergerate_pop-pp.decayrate_pop)*vv.pop_growth_region[tt-1, rr]+pp.decayrate_pop*sum(vv.pop_growth_region[tt-1, :] .* pp.popweights_region)
                end
            end

            for cc in dd.country
                rr = vv.country2region[cc]
                if rr == 0
                    growth = 0
                else
                    growth = vv.pop_growth_region[tt, rr]
                end

                vv.pop[tt, cc] = vv.pop[tt-1, cc] * (1 + growth)
            end
        end

        vv.utility_sum[tt] = 0
        vv.world_population[tt] = 0

        for cc in dd.country

            if pp.conspc[tt, cc] == 0

               vv.utility[tt, cc] = 0

            else

               vv.utility[tt, cc] = (1 / (1 - pp.EMUC) * (max(pp.conspc[tt, cc], 1) * pp.lossfactor[tt, cc]) ^ (1 - pp.EMUC)) * vv.pop[tt, cc]
            end

            vv.utility_sum[tt] += vv.utility[tt, cc]
            vv.world_population[tt] += vv.pop[tt, cc]

        end

        vv.world_disc_utility[tt] = vv.utility_sum[tt] * (1 + pp.PRTP) ^ - (gettime(tt) - 2020) # possibly scope to improve how time is coded in exponent
        vv.equiv_conspc[tt] = (vv.world_disc_utility[tt] * (1 - pp.EMUC) / vv.world_population[tt]) ^ (1 / (1 - pp.EMUC))
    end
end

function addUtility(model, ssp)

    params = CSV.read("../data/utilityparams.csv", DataFrame)

    utility = add_comp!(model, Utility, first=2010)

    utility[:EMUC] = params.Value[params.Parameter .== "EMUC"][1]
    utility[:PRTP] = params.Value[params.Parameter .== "PRTP"][1]

    utility[:ssp] = ssp

    sspextend = CSV.read("../data/sspextend-pop.csv", DataFrame)
    utility[:convergerate_pop] = sspextend[sspextend.SSP .== ssp, "Convergence rate"][1]
    utility[:decayrate_pop] = sspextend[sspextend.SSP .== ssp, "Decay rate"][1]
    utility[:popweights_region] = [getpopweight_ssp(region, ssp) for region in dim_keys(model, :region)]

    utility[:lossfactor] = reshape(ones(dim_count(model, :time) * dim_count(model, :country)),
                                   dim_count(model, :time), dim_count(model, :country))

    utility

end
