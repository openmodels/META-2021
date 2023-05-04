@defcomp TemperatureConverter begin
    # Variables
    T_AT = Variable(index=[time]) # after bias-correction
    biascorrection = Variable()

    # Parameters
    T = Parameter(index=[time]) # raw from MimiFAIR

    function run_timestep(pp, vv, dd, tt)
        if gettime(tt) == 2005
            # Fill in all values
            vv.biascorrection = 0.69 - mean([pp.T[tt-dt] for dt in 0:19]) # relative to 1850-1900
            ss = 1
            while TimestepIndex(ss) <= tt
                vv.T_AT[TimestepIndex(ss)] = pp.T[TimestepIndex(ss)] + vv.biascorrection
                ss += 1
            end
        elseif gettime(tt) > 2005
            vv.T_AT[tt] = pp.T[tt] + vv.biascorrection
        end
    end
end

function addTemperatureConverter(model; before=nothing, after=nothing)
    add_comp!(model, TemperatureConverter, before=before, after=after)
end
