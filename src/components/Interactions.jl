include("../lib/interacts.jl")

@defcomp Interactions begin
    # Variables
    p_GIS = Variable(index=[time])
    p_NINO = Variable(index=[time])
    p_AIS = Variable(index=[time]) # NEW

    f_AMOC = Variable(index=[time])
    f_GIS = Variable(index=[time])
    f_WAIS = Variable(index=[time])
    f_AMAZ = Variable(index=[time])
    f_NINO = Variable(index=[time])
    f_AIS = Variable(index=[time]) # NEW

    # Parameters
    I_AMOC = Parameter(index=[time])
    VGIS = Parameter(index=[time], unit="percent")
    I_WAIS = Parameter(index=[time])
    I_AMAZ = Parameter(index=[time])
    mNINO3pt4 = Parameter(index=[time], unit="MSLP")
    totalSLR_Ross = Parameter(index=[time], unit="m") # NEW
    totalSLR_Amundsen = Parameter(index=[time], unit="m") # NEW
    totalSLR_Weddell = Parameter(index=[time], unit="m") # NEW
    totalSLR_Peninsula = Parameter(index=[time], unit="m") # NEW

    gis2amoc = Parameter()
    wais2amoc = Parameter()
    amaz2amoc = Parameter()
    nino2amoc = Parameter()
    amoc2gis = Parameter()
    wais2gis = Parameter()
    amaz2gis = Parameter()
    nino2gis = Parameter()
    amoc2wais = Parameter()
    gis2wais = Parameter()
    amaz2wais = Parameter()
    nino2wais = Parameter()
    amoc2amaz = Parameter()
    gis2amaz = Parameter()
    wais2amaz = Parameter()
    nino2amaz = Parameter()
    amoc2nino = Parameter()
    gis2nino = Parameter()
    wais2nino = Parameter()
    amaz2nino = Parameter()

    function run_timestep(pp, vv, dd, tt)
        if is_first(tt)
            vv.f_AMOC[tt] = 1
            vv.f_GIS[tt] = 1
            vv.f_WAIS[tt] = 1
            vv.f_AMAZ[tt] = 1
            vv.f_NINO[tt] = 1
            vv.f_AIS[tt] = 1
        else
            vv.p_GIS[tt-1] = 1 - pp.VGIS[tt-1]
            vv.p_NINO[tt-1] = max(min(10000 * (1 - pp.mNINO3pt4[tt-1] / pp.mNINO3pt4[TimestepIndex(1)]), 1), 0)
            vv.p_AIS[tt-1] = 1/5 * ( pp.totalSLR_Ross[tt-1] + pp.totalSLR_Amundsen[tt-1] + pp.totalSLR_Weddell[tt-1] + pp.totalSLR_Peninsula[tt-1] ) # NEW

            ## Template: (pp.amoc2foo^pp.I_AMOC[tt-1]) * (pp.gis2foo^vv.p_GIS[tt-1]) * (pp.wais2foo^pp.I_WAIS[tt-1]) * (pp.amaz2foo^pp.I_AMAZ[tt-1]) * (pp.nino2foo^vv.p_NINO[tt-1])
            ## Note WAIS interaction can be via old WAIS module or new AIS-based estimate of interaction
            vv.f_AMOC[tt] = (pp.gis2amoc^vv.p_GIS[tt-1]) * (pp.wais2amoc^pp.I_WAIS[tt-1]) * (pp.amaz2amoc^pp.I_AMAZ[tt-1]) * (pp.nino2amoc^vv.p_NINO[tt-1]) * (pp.wais2amoc^vv.p_AIS[tt-1])
            vv.f_GIS[tt] = (pp.amoc2gis^pp.I_AMOC[tt-1]) * (pp.wais2gis^pp.I_WAIS[tt-1]) * (pp.amaz2gis^pp.I_AMAZ[tt-1]) * (pp.nino2gis^vv.p_NINO[tt-1]) * (pp.wais2gis^vv.p_AIS[tt-1])
            vv.f_WAIS[tt] = (pp.amoc2wais^pp.I_AMOC[tt-1]) * (pp.gis2wais^vv.p_GIS[tt-1]) * (pp.amaz2wais^pp.I_AMAZ[tt-1]) * (pp.nino2wais^vv.p_NINO[tt-1])
            vv.f_AMAZ[tt] = (pp.amoc2amaz^pp.I_AMOC[tt-1]) * (pp.gis2amaz^vv.p_GIS[tt-1]) * (pp.wais2amaz^pp.I_WAIS[tt-1]) * (pp.nino2amaz^vv.p_NINO[tt-1]) * (pp.wais2amaz^vv.p_AIS[tt-1])
            vv.f_NINO[tt] = (pp.amoc2nino^pp.I_AMOC[tt-1]) * (pp.gis2nino^vv.p_GIS[tt-1]) * (pp.wais2nino^pp.I_WAIS[tt-1]) * (pp.amaz2nino^pp.I_AMAZ[tt-1]) * (pp.wais2nino^vv.p_AIS[tt-1])
            vv.f_AIS[tt] = (pp.amoc2wais^pp.I_AMOC[tt-1]) * (pp.gis2wais^vv.p_GIS[tt-1]) * (pp.amaz2wais^pp.I_AMAZ[tt-1]) * (pp.nino2wais^vv.p_NINO[tt-1]) # NEW # note AIS probability factors same as WAIS
        end
    end
end

function addInteractions(model; before=nothing, after=nothing)
    interact = add_comp!(model, Interactions, first=2010, before=before, after=after)

    allinteractrates((symbol, ratemu, ratese) -> interact[symbol] = ratemu)
    interact[:I_AMOC] = zeros(dim_count(model, :time))
    interact[:VGIS] = ones(dim_count(model, :time))
    interact[:I_WAIS] = zeros(dim_count(model, :time))
    interact[:I_AMAZ] = zeros(dim_count(model, :time))
    interact[:mNINO3pt4] = ones(dim_count(model, :time))
    interact[:totalSLR_Ross] = zeros(dim_count(model, :time)) # NEW
    interact[:totalSLR_Amundsen] = zeros(dim_count(model, :time)) # NEW
    interact[:totalSLR_Weddell] = zeros(dim_count(model, :time)) # NEW
    interact[:totalSLR_Peninsula] = zeros(dim_count(model, :time)) # NEW

    interact
end
