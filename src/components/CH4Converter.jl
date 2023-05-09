@defcomp CH4Converter begin
    # Variables
    E_ch4 = Variable(index=[time])   # Annual methane emissions (TgCH₄ yr⁻¹)

    # Parameters #ALL READ IN FROM RCP SCENARIOS AND OTHER TP MODULES
    ch4_2009 = Parameter(unit="MtCH4", default=364.2003634)
    ch4_rcp = Parameter(index=[time], unit="MtCH4") # Named based on James's variable from RCP.jl
    ch4_pcf = Parameter(index=[time], unit="MtCH4")
    ch4_omh = Parameter(index=[time], unit="MtCH4")
    ch4_extra = Parameter(index=[time], unit="MtCH4")

    function run_timestep(pp, vv, dd, tt)
        if is_first(tt)
            vv.E_ch4[tt] = pp.ch4_2009
        else
            vv.E_ch4[tt] = pp.ch4_rcp[tt-1] + pp.ch4_pcf[tt-1] + pp.ch4_omh[tt-1] + pp.ch4_extra[tt-1]
        end
    end
end

function addCH4Converter(model; before=nothing, after=nothing)
    ch4converter = add_comp!(model, CH4Converter, first=2010, before=before, after=after)

    #Prefill CH4 TPs with zeros until we have the components defined (change later once we link to PCF and OMH modules)
    ch4converter[:ch4_pcf] = zeros(dim_count(model, :time))
    ch4converter[:ch4_omh] = zeros(dim_count(model, :time))
    ch4converter[:ch4_extra] = zeros(dim_count(model, :time))

    ch4converter
end
