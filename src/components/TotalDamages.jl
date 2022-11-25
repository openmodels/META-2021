@defcomp TotalDamages begin

        # Per capita costs to aggregate
        postdamage_consumption_percap_percountry  = Parameter(index = [time, country], unit = "2010 USD PPP") # Consumption per cap per country post SLR post temp damages
        baseline_consumption_percap_percountry    = Parameter(index = [time, country], unit = "2010 USD PPP") # Counterfactual consumption per cap per country from SSPs
        market_damages_percap_peryear             = Variable(index = [time, country], unit = "2010 USD PPP") # Market damages per capita per year

        # Other parameters
        population = Parameter(index = [time, country], unit = "million person")
        #period_length = Parameter(index = [time], unit="year") #Can use 2050 (30 year length) or 2100 (80 year length). However, can I just multiply yearly damages by amount of years? Shouldn't this be a sum, since yearly damages are different for each year?

        # Undiscounted costs to calculate
        total_damages_percap_peryear    = Variable(index = [time, country], unit = "2010 USD PPP")  # Undiscounted total damages per capita per year   
        total_damages_peryear           = Variable(index = [time, country], unit = "Million 2010 USD PPP")  # Undiscounted total damages per year per country
        total_damages_global_peryear    = Variable(index = [time], unit = "Million 2010 USD PPP")  # Undiscounted total damages per year 
        total_damages_global_cumulative = Variable(index = [time], unit = "Million 2010 USD PPP") # Undiscounted total damages through 2200

        function run_timestep(pp, vv, dd, tt)
                #Define market damages per cap per year
                for cc in dd.country
                        vv.market_damages_percap_peryear[tt,cc] = pp.baseline_consumption_percap_percountry[tt,cc] - pp.postdamage_consumption_percap_percountry[tt,cc]
                end

                # Undiscounted damages
                vv.total_damages_percap_peryear[tt, :] = pp.market_damages_percap_peryear[t, :] # This is not really needed, but perhaps good in case we add new damages channels at some point.
                vv.total_damages_peryear[tt, :] = vv.total_damages_percap_peryear[tt, :] .* pp.population[tt, :]
                vv.total_damages_global_peryear[tt] = 0

                for cc in dd.country
                        vv.total_damages_global_peryear[tt] += vv.total_damages_peryear[tt, cc]
                end

                if is_first(tt)
                        vv.total_damages_global_cumulative[tt] = vv.total_damages_global_peryear[tt]
                else
                        vv.total_damages_global_cumulative[tt] = vv.total_damages_global_cumulative[tt-1] + vv.total_damages_global_peryear[tt]
                end
        end

end


function addTotalDamages(model)

        params = CSV.read("../data/utilityparams.csv", DataFrame)
    
        damages = add_comp!(model, TotalDamages)
    
        damages
end
