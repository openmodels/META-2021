@defcomp TemperatureConverter begin
    # Variables
    T_AT = Variable(index=[time]) # after bias-correction
    T_AT_tminus100 = Variable(index=[time]) # T_AT 100 years before
    biascorrection = Variable()

    # Parameters
    T = Parameter(index=[time]) # raw from MimiFAIR

    function run_timestep(pp, vv, dd, tt)
        if gettime(tt) == 2010
            # Fill in all values
            vv.biascorrection = 0.854 - pp.T[tt]
            ss = 1
            while TimestepIndex(ss) <= tt
                vv.T_AT[TimestepIndex(ss)] = pp.T[TimestepIndex(ss)] + vv.biascorrection
                ss += 1
            end
            vv.T_AT_tminus100[tt] = pp.T_AT[tt - 100]
        elseif gettime(tt) > 2010
            vv.T_AT[tt] = pp.T[tt] + vv.biascorrection
            vv.T_AT_tminus100[tt] = pp.T_AT[tt - 100]
        else
            vv.T_AT_tminus100[tt] = 0
        end
    end
end

function addTemperatureConverter(model; before=nothing, after=nothing)
    add_comp!(model, TemperatureConverter, before=before, after=after)
end
