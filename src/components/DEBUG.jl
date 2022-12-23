@defcomp DEBUG begin

    # Declare variable
    utility_counterfactual              = Variable(index =[time, region])# Utility per country of counterfactual no climate damages SSP consumption

    function run_timestep(pp, vv, dd, tt)
        # Compute variable
        for cc in dd.country
            vv.utility_counterfactual[tt, cc] = 0
        end
    end
end

function addDEBUG(model)

        debug = add_comp!(model, DEBUG)

        debug
end