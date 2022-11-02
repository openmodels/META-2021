@defcomp WAISmodel begin
    # Variables
    p_WAIS = Variable(index=[time])
    I_WAIS = Variable(index=[time])
    contemp_SLR_WAIS = Variable(index=[time], unit="m")
    SLR_WAIS = Variable(index=[time], unit="m")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC") # Link manually to other input later
    uniforms = Parameter(index=[time])

    b_WAIS = Parameter(default=0.0043) # Hazard rate
    waisrate = Parameter(unit="m/year", default=0.0033) # Melting rate (in m/year)
    #= This is the mean melting rate used in Diaz and Keller (2016). META2021 can also use their 2.5%ile/97.5%ile,
    whereas the original paper mentions that the rate is drawn from a log-normal distribution with mean=0.0033 and
    standard deviation 0.00165. We could implement this here, but so far I have not done it to reproduce META2021
    before making changes.=#

    f_WAIS = Parameter(index=[time])

    function run_timestep(pp, vv, dd, tt)
        if is_first(tt)
            vv.p_WAIS[tt] = 0 # Likelihood of WAIS disintegration being triggered
            vv.I_WAIS[tt] = 0 # Indicator function for WAIS disintegration is triggered
            vv.SLR_WAIS[tt] = 0 # Sea-level rise from WAIS melting (in m)

        else
            vv.p_WAIS[tt] = min(pp.f_WAIS[tt] * (pp.b_WAIS*pp.T_AT[tt]-0.6)^2, 1)

            if vv.I_WAIS[tt-1] != 1
                vv.I_WAIS[tt] = pp.uniforms[tt] < vv.p_WAIS[tt]
                #pp.I_WAIS[tt] = rand(Binomial(1,vv.p_WAIS[tt]))

            else
                vv.I_WAIS[tt] = 1

            end

            if vv.I_WAIS[tt] != 1
                vv.contemp_SLR_WAIS[tt] = 0
                vv.SLR_WAIS[tt] = vv.SLR_WAIS[tt-1]
            else
                vv.contemp_SLR_WAIS[tt] = pp.waisrate
                vv.SLR_WAIS[tt] = vv.SLR_WAIS[tt-1] + vv.contemp_SLR_WAIS[tt]

            end

        end
    end
end

function addWAISmodel(model, default; before=nothing, after=nothing)

    #if WAIScalib == "Distribution"
    #    error("Distribution WAIS model not implemented")
    #end

    waismodel = add_comp!(model, WAISmodel, before=before, after=after)
    waismodel[:f_WAIS] = ones(dim_count(model, :time))
    waismodel
end
