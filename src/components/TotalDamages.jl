using Mimi
@defcomp TotalDamages begin

        # Per capita values for aggregation
        postdamage_consumption_percap_percountry  = Parameter(index = [time, country], unit = "2010 USD PPP") # Consumption per cap per country post SLR post temp damages
        baseline_consumption_percap_percountry    = Parameter(index = [time, country], unit = "2010 USD PPP") # Counterfactual consumption per cap per country from SSPs
        market_damages_percap_peryear             = Variable(index = [time, country], unit = "2010 USD PPP") # Market damages per capita per year
        market_damages_percap_peryear_percent     = Variable(index = [time, country], unit = "Per cent") # % of consumption per cap lost due to market damages from climate change
        global_conspc_counterfactual              = Variable(index = [time], unit = "2010 USD PPP") #Global consumption per capita before climate damages
        total_damages_equiv_conspc_equity         = Variable(index = [time]) # Per period consumption-equivalent for equity-weighted climate damages

        # Other parameters and variables
        population                                = Parameter(index = [time, country], unit = "inhabitants")
        population_global                         = Variable(index = [time], unit = "inhabitants")
        population_weights                        = Variable(index = [time, country], unit = "Per cent")
        EMUC                                      = Parameter() # elasticity of marginal utility of consumption
        lossfactor                                = Parameter(index=[time, country]) # non-market loss factor for MERGE add-on
        lossfactor_global                         = Variable(index = [time]) # Global loss-factor, from population-weighted country loss-factors

        # Undiscounted items to calculate
        total_damages_percap_peryear              = Variable(index = [time, country], unit = "2010 USD PPP")  # Undiscounted total damages per capita per year per country
        total_damages_percap_peryear_percent      = Variable(index = [time, country])  # Undiscounted total damages per capita per year per country in % of baseline consumption
        total_damages_peryear                     = Variable(index = [time, country], unit = "2010 USD PPP")  # Undiscounted total damages per year per country
        total_damages_cumulative                  = Variable(index = [time, country], unit = "2010 USD PPP") # Country-level undiscounted total damages through 2200
        total_damages_global_peryear              = Variable(index = [time], unit = "2010 USD PPP")  # Undiscounted total damages per year
        total_damages_global_peryear_percent      = Variable(index = [time]) # Undiscounted total damages per year in % of baseline consumption
        total_damages_global_cumulative           = Variable(index = [time], unit = "2010 USD PPP") # Undiscounted total damages through 2200
        utility_equivalent_change                 = Variable(index = [time, country]) # Welfare-equivalent change from climate damages for change in consumption per capita due to climate damages per country 
        utility_equivalent_change_global          = Variable(index = [time]) # Welfare-equivalent change from climate damages population-weighted (accounts for equity weighting through utility function)
        total_damages_equiv_conspc_equity         = Variable(index = [time]) # Per period consumption-equivalent for equity-weighted climate damages
        

        function run_timestep(pp, vv, dd, tt)
                #Define market damages per cap per year
                for cc in dd.country
                        vv.market_damages_percap_peryear[tt,cc] = pp.baseline_consumption_percap_percountry[tt,cc] - pp.postdamage_consumption_percap_percountry[tt,cc]
                        vv.market_damages_percap_peryear_percent[tt,cc] = vv.market_damages_percap_peryear[tt,cc]/pp.baseline_consumption_percap_percountry[tt,cc]*100
                        #Utility equivalent change per country per person
                        vv.utility_equivalent_change[tt, cc] = ((1 / (1 - pp.EMUC)*(max(pp.baseline_consumption_percap_percountry[tt, cc],1)*pp.lossfactor[tt, cc])^(1 - pp.EMUC))) - ((1 / (1 - pp.EMUC)*(max(pp.postdamage_consumption_percap_percountry[tt,cc],1)*pp.lossfactor[tt, cc])^(1 - pp.EMUC)))
                end
        
                # Undiscounted damages
                vv.total_damages_percap_peryear[tt, :] = vv.market_damages_percap_peryear[tt, :] # This is not really needed, but perhaps good in case we add new damages channels at some point.
                vv.total_damages_percap_peryear_percent[tt, :] = vv.market_damages_percap_peryear_percent[tt, :]
              
                vv.total_damages_global_peryear[tt] = 0
                vv.population_global[tt] = 0
                vv.total_damages_global_peryear_percent[tt] = 0
                vv.utility_equivalent_change_global[tt] = 0
                vv.global_conspc_counterfactual[tt] = 0

                for cc in dd.country
                        if isnan(pp.population[tt, cc])
                                pp.population[tt, cc] = 0
                        end
                        vv.total_damages_peryear[tt, cc] = vv.total_damages_percap_peryear[tt, cc] * pp.population[tt, cc]
                        vv.total_damages_global_peryear[tt] += vv.total_damages_peryear[tt, cc]
                        vv.population_global[tt] += pp.population[tt, cc]
                        #Sanitize for missing countries
                        if isnan(vv.total_damages_percap_peryear_percent[tt, cc])
                                vv.total_damages_percap_peryear_percent[tt, cc] = 0
                        end
                end

                for cc in dd.country
                        
                        #The variable below computes population-weighted percentage-change per year
                        vv.population_weights[tt, cc] = pp.population[tt, cc] / vv.population_global[tt]
                        vv.total_damages_global_peryear_percent[tt] += vv.total_damages_percap_peryear_percent[tt, cc] * vv.population_weights[tt, cc]
                        vv.utility_equivalent_change_global[tt] += vv.utility_equivalent_change[tt,cc] * vv.population_weights[tt, cc]
                        #Global consumption per capita pop-weighted
                        vv.global_conspc_counterfactual[tt] += pp.baseline_consumption_percap_percountry[tt,cc]*vv.population_weights[tt, cc]
                
                end
                
                #Transform equity-weighted number in consumption equivalent in % change
                vv.lossfactor_global[tt] = 0
                for cc in dd.country
                        vv.lossfactor_global[tt] += pp.lossfactor[tt, cc]*(pp.population[tt, cc]/vv.population_global[tt])
                end
                vv.total_damages_equiv_conspc_equity[tt] = (vv.global_conspc_counterfactual[tt]-((vv.utility_equivalent_change_global[tt]*vv.lossfactor_global[tt]*(1 - pp.EMUC)) ^ (1 / (1 - pp.EMUC))))/vv.global_conspc_counterfactual[tt]
                
                if is_first(tt)
                        vv.total_damages_global_cumulative[tt] = vv.total_damages_global_peryear[tt]
                        for cc in dd.country
                                vv.total_damages_cumulative[tt, cc] = vv.total_damages_peryear[tt, cc] # Initializes cumulative damages per country
                        end
                else
                        vv.total_damages_global_cumulative[tt] = vv.total_damages_global_cumulative[tt-1] + vv.total_damages_global_peryear[tt]
                        for cc in dd.country
                                vv.total_damages_cumulative[tt, cc] = vv.total_damages_cumulative[tt-1, cc] + vv.total_damages_peryear[tt,cc] # Need to sum country-level total damages for country-level results
                        end
                end
                
        end

end


function addTotalDamages(model)

        params = CSV.read("../data/utilityparams.csv", DataFrame)

        damages = add_comp!(model, TotalDamages, first=2010)
        damages[:EMUC] = params.Value[params.Parameter .== "EMUC"][1]
        damages[:lossfactor] = reshape(ones(dim_count(model, :time) * dim_count(model, :country)),
        dim_count(model, :time), dim_count(model, :country))
        
        damages
end
