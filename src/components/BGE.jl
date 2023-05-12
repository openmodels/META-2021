@defcomp BGE begin

    # Balanced growth equivalent inputs
    world_utility                           = Variable(index =[time])# World utility
    world_utility_counterfactual            = Variable(index =[time])# World utility of counterfactual no climate damages SSP consumption
    world_disc_utility                      = Variable(index =[time]) # Calculate myself similar to utility module, but with variable cutoff year (2050).
    world_disc_utility_counterfactual       = Variable(index =[time]) # Calculates utility for SSP scenario consumption without climate damages
    sum_world_disc_utility                  = Variable() # Sum discounted world utility at present. Scalar
    sum_world_disc_utility_counterfactual   = Variable() # Sum discounted counterfactual world utility at present. Scalar

    utility                                 = Parameter(index =[time, country])# Utility per country
    utility_counterfactual                  = Variable(index =[time, country])# Utility per country of counterfactual no climate damages SSP consumption
    disc_utility                            = Variable(index =[time, country])
    disc_utility_counterfactual             = Variable(index =[time, country])
    sum_disc_utility                        = Variable(index =[country]) # Discounted utility at present. One-dimensional array
    sum_disc_utility_counterfactual         = Variable(index =[country]) # Discounted utility at present. One-dimensional array

    pop                                     = Parameter(index=[time, country], unit="inhabitants") # national population

    # Utility parameters
    EMUC                                    = Parameter() # elasticity of marginal utility of consumption
    PRTP                                    = Parameter(unit ="pct") # pure rate of time preference
    lossfactor                              = Parameter(index =[time, country]) # non-market loss factor


    # Consumption variables
    baseline_consumption_percap_percountry = Parameter(index = [time, country], unit = "2010 USD PPP") # Counterfactual consumption per cap per country from SSPs

    # Balanced-growth equivalent outputs
    world_bge                           = Variable()
    bge                                 = Variable(index =[country])

    function run_timestep(pp, vv, dd, tt)
        # Set up empty world utility var buckets
        vv.world_utility_counterfactual[tt] = 0
        vv.world_utility[tt] = 0

        #Calculate counterfactual utility of SSP no climate damages-scenario
        for cc in dd.country

            if pp.baseline_consumption_percap_percountry[tt, cc] == 0

                vv.utility_counterfactual[tt, cc] = 0

            else

                #Cut damage calculations in 2050 for methane paper
                if tt <= TimestepValue(2050)

                    vv.utility_counterfactual[tt, cc] = (1 / (1 - pp.EMUC) * (max(pp.baseline_consumption_percap_percountry[tt, cc], 1) * pp.lossfactor[tt, cc]) ^ (1 - pp.EMUC)) * pp.pop[tt, cc]

                else
                    vv.utility_counterfactual[tt, cc] = 0
                    pp.utility[tt, cc] = 0

                end
            end

            # Calculate world utility
            vv.world_utility_counterfactual[tt] += vv.utility_counterfactual[tt, cc]
            vv.world_utility[tt] += pp.utility[tt, cc]
        end

        # Discount utility
        vv.world_disc_utility_counterfactual[tt] = vv.world_utility_counterfactual[tt] * (1 + pp.PRTP) ^ - (gettime(tt) - 2020)
        vv.world_disc_utility[tt] = vv.world_utility[tt] * (1 + pp.PRTP) ^ - (gettime(tt) - 2020)

        for cc in dd.country
            vv.disc_utility_counterfactual[tt, cc] = vv.utility_counterfactual[tt, cc] * (1 + pp.PRTP) ^ - (gettime(tt) - 2020)
            vv.disc_utility[tt, cc] = pp.utility[tt, cc] * (1 + pp.PRTP) ^ - (gettime(tt) - 2020)
        end

                #= Marginal welfare change method
        vv.world_damages_marginalmethod[tt] = (world_disc_utility_counterfactual[tt] - world_disc_utility[tt])#(1/marginal utility of cosumption in 2010^-EMUC) WHERE TO FIND CONSUMPTION IN 2010? AND CAN I USE CONSUMPTION PER CAP?
        vv.damages_marginalmethod[tt, cc] = (disc_utility_counterfactual[tt, cc] - disc_utility[tt, cc])#(1/marginal utility of cosumption in 2010^-EMUC) WHERE TO FIND CONSUMPTION IN 2010? AND CAN I USE CONSUMPTION PER CAP?
        =#


    end
end

function addBGE(model)

    params = CSV.read("../data/utilityparams.csv", DataFrame)

    bge_comp = add_comp!(model, BGE, first=2010)

    #Grabbing these from Utility.jl via MimiMETA.jl somehow did not work
    bge_comp[:EMUC] = params.Value[params.Parameter .== "EMUC"][1]
    bge_comp[:PRTP] = params.Value[params.Parameter .== "PRTP"][1]
    bge_comp[:lossfactor] = reshape(ones(dim_count(model, :time) * dim_count(model, :country)),
    dim_count(model, :time), dim_count(model, :country))

    bge_comp
end
