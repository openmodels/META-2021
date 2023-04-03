@defcomp CO2Converter begin
    # Variables
    E_co2 = Variable(index=[time])   # Annual carbon dioxide emissions (GtC yr⁻¹).

    # Parameters
    co2_2009 = Parameter(unit="GtCO2", default=35.70037753)
    co2_rcp = Parameter(index=[time], unit="GtCO2")
    co2_pcf = Parameter(index=[time], unit="GtCO2")
    co2_amazon = Parameter(index=[time], unit="GtCO2")
    co2_extra = Parameter(index=[time], unit="GtCO2")

    function run_timestep(pp, vv, dd, tt)
        if is_first(tt)
            vv.E_co2[tt] = pp.co2_2009 / 3.67
        else
            vv.E_co2[tt] = (pp.co2_rcp[tt-1] / 1000 + pp.co2_pcf[tt-1] + pp.co2_amazon[tt-1] + pp.co2_extra[tt]) / 3.67
        end
    end
end

function addCO2Converter(model; before=nothing, after=nothing)
    co2converter = add_comp!(model, CO2Converter, first=2010, before=before, after=after)

    co2converter[:co2_pcf] = zeros(dim_count(model, :time))
    co2converter[:co2_amazon] = zeros(dim_count(model, :time))
    co2converter[:co2_extra] = zeros(dim_count(model, :time))

    co2converter
end
