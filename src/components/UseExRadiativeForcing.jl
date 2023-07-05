@defcomp UseExRadiativeForcing begin
    # Variables
    total_RF = Variable(index=[time])

    # Parameters
    co2_RF = Parameter(index=[time])
    ch4_RF = Parameter(index=[time])
    ch4_o3_RF = Parameter(index=[time])
    ch4_h2o_RF = Parameter(index=[time])

    F_EX = Parameter(index=[time])

    total_RF_FAIR = Parameter(index=[time])

    function run_timestep(pp, vv, dd, tt)
        if gettime(tt) >= 2010
            vv.total_RF[tt] = pp.co2_RF[tt] + pp.ch4_RF[tt] +
                pp.ch4_o3_RF[tt] + pp.ch4_h2o_RF[tt] + pp.F_EX[tt]
        else
            vv.total_RF[tt] = pp.total_RF_FAIR[tt]
        end
    end
end

function addUseExRadiativeForcing(model; before=nothing, after=nothing)
    add_comp!(model, UseExRadiativeForcing, before=before, after=after)
end
