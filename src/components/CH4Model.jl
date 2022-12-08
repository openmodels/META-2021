@defcomp CH4Model begin
    # Variables
    CH4_concentration = Variable(index=[time], unit="ppb")
    I_CH4_N2O = Variable(index=[time])
    F_CH4 = Variable(index=[time], unit="W/m^2")
    N2O_initial_concentration = Variable(index=[time], unit="ppb")
    fMN_2010 = Variable(index=[time])

    # Parameters #ALL READ IN FROM RCP SCENARIOS AND OTHER TP MODULES
    ch4_rcp = Parameter(index=[time], unit="MtCH4") # Named based on James's variable from RCP.jl
    ch4_pcf = Parameter(index=[time], unit="MtCH4")
    ch4_omh = Parameter(index=[time], unit="MtCH4")
    ch4_extra = Parameter(index=[time], unit="MtCH4")
    ch4_conc_rcp = Parameter(index=[time], unit="ppb")
    n2o_conc_rcp = Parameter(index=[time], unit="ppb")

    ch4_alpha = Parameter(default=0.036)
    fMN_parameter1 = Parameter(default=0.0000201)
    fMN_parameter2 = Parameter(default=0.00000000000000531)
    decay_rate = Parameter(default=1/11.8) # Upgraded to 11.8 years of atmospheric residence time to reflect AR6. Previous AR5 value used in Dietz et al. 2021 PNAS is 12.4.
    conversion_ppb_Mt = Parameter(default=2.78)
    CH4_conc_preindustrial = Parameter(unit="ppb", default=722)


function run_timestep(pp, vv, dd, tt)
        if is_first(tt)
            vv.CH4_concentration[tt] = pp.ch4_conc_rcp[TimestepIndex(1)]
            vv.N2O_initial_concentration[tt] = pp.n2o_conc_rcp[TimestepIndex(1)]
            vv.fMN_2010[tt] = 0.47*log(1+pp.fMN_parameter1*(pp.CH4_conc_preindustrial*vv.N2O_initial_concentration[TimestepIndex(1)])^0.75+pp.fMN_parameter2*pp.CH4_conc_preindustrial*(pp.CH4_conc_preindustrial*vv.N2O_initial_concentration[TimestepIndex(1)])^1.52) #James suggested to put this calculation into an init() function that precedes run_timestep(), but I need the latter to grab the first element of the N2O concentration array.

        else

        vv.CH4_concentration[tt] = pp.CH4_conc_preindustrial +(vv.CH4_concentration[tt-1]-pp.CH4_conc_preindustrial)*(1-pp.decay_rate)+(pp.ch4_rcp[tt-1]+pp.ch4_pcf[tt-1]+pp.ch4_omh[tt-1]+pp.ch4_extra[tt-1])/pp.conversion_ppb_Mt

            if vv.CH4_concentration[tt] < 0
                vv.CH4_concentration[tt] = 0
            end
        end

        vv.I_CH4_N2O[tt] = 0.47*log(1+pp.fMN_parameter1*(vv.CH4_concentration[tt]*vv.N2O_initial_concentration[TimestepIndex(1)])^(0.75)+pp.fMN_parameter2*vv.CH4_concentration[tt]*(vv.CH4_concentration[tt]*vv.N2O_initial_concentration[TimestepIndex(1)])^(1.52))
        vv.F_CH4[tt] = pp.ch4_alpha*(sqrt(vv.CH4_concentration[tt])-sqrt(pp.CH4_conc_preindustrial))-(vv.I_CH4_N2O[tt]-vv.fMN_2010[TimestepIndex(1)])
    end
end

function addCH4Model(model, ch4calib)
    ch4model = add_comp!(model, CH4Model)

    if ch4calib == "default"
        ch4model[:ch4_alpha] = 0.036

    elseif ch4calib == "low"
        ch4model[:ch4_alpha] = 0.0319967

    elseif ch4calib == "high"
        ch4model[:ch4_alpha] = 0.0400033

    elseif ch4calib == "Distribution"
        error("Distribution CH4 cycle not implemented")
    end

    #Prefill CH4 TPs with zeros until we have the components defined (change later once we link to PCF and OMH modules)
    ch4model[:ch4_pcf] = zeros(dim_count(model, :time))
    ch4model[:ch4_omh] = zeros(dim_count(model, :time))
    ch4model[:ch4_extra] = zeros(dim_count(model, :time))

    ch4model
end
