@defcomp SLRModel begin
    # Variables
    SLR_therm = Variable(index=[time], unit="m")
    SLR = Variable(index=[time], unit="m")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC")
    SLR_GIS = Parameter(index=[time], unit="m")
    SLR_WAIS = Parameter(index=[time], unit="m")
    SLR_AIS = Parameter(index=[time], unit="m")

    r_TE = Parameter(default=0.000779) #SLR from thermal expansion
    r_GSIC = Parameter(default=0.0008164) #SLR from glaciers and small ice caps (Excel META: =0.0314/10*0.26)

    SLR_initial = Parameter(default=0.04, unit="m")

    function run_timestep(pp, vv, dd, tt)
        if is_first(tt)
            vv.SLR_therm[tt] = pp.SLR_initial
            vv.SLR[tt] = pp.SLR_initial
        else
            # Calculate thermal SLR contribution
            vv.SLR_therm[tt] = (pp.r_TE + pp.r_GSIC)*pp.T_AT[tt] + vv.SLR_therm[tt-1]

            # Calculate total SLR
            vv.SLR[tt] = vv.SLR_therm[tt] + pp.SLR_GIS[tt] + pp.SLR_WAIS[tt] + pp.SLR_AIS[tt]
        end
    end
end

function addSLR(model)
    slr = add_comp!(model, SLRModel, first=2010)

    slr[:SLR_GIS] = zeros(dim_count(model, :time))
    slr[:SLR_WAIS] = zeros(dim_count(model, :time))
    slr[:SLR_AIS] = zeros(dim_count(model, :time))

    slr
end
