@defcomp GISModel begin

    # Variables

    VGIS = Variable(index=[time], unit="percent") # volume of the GIS relative to initial volume
    Diff_VGIS = Variable(index=[time], unit="percent") # difference in VGIS between periods
    T_GISstar = Variable(index=[time], unit="degC") # temperature increase relative to initial temperature associated with a particular degree of GIS melting in equilibrium
    TD = Variable(index=[time], unit="degC") # difference between T_AT and T_GISstar
    SLR_GIS = Variable(index=[time], unit="metres") # sea-level rise from GIS melting

    # Parameters

    T_AT = Parameter(index=[time], unit="degC") # Atmospheric temperature relative to pre-industrial, output of temperature model

    meltmult = Parameter(default=1) # Multiplier times melt rate
    volzero = Parameter(default=0.999, unit="percent") # Initial volume of GIS (%)
    avoldot = Parameter() # Equilibrium volume temperature equation
    tmaxa = Parameter(unit="degC") # Global temperature at which minimum volume
    slrgis = Parameter(default=7, unit="metres") # Sea level rise from GIS melt (metres)
    expvol = Parameter(default=0.2) # Exponent on voldot
    exptstar = Parameter() # Exponent on Tstar

    f_GIS = Parameter(index=[time])

    function run_timestep(pp, vv, dd, tt)

        if is_first(tt)

           vv.VGIS[tt] = pp.volzero
           vv.T_GISstar[tt] = pp.tmaxa * (1 - vv.VGIS[tt]) ^ pp.exptstar
           vv.TD[tt] = pp.T_AT[tt] - vv.T_GISstar[tt] + 0.00001

        else

           vv.Diff_VGIS[tt] = pp.avoldot * sign(vv.TD[tt-1]) * vv.TD[tt-1] ^ 2 * vv.VGIS[tt-1] ^ pp.expvol
           vv.Diff_VGIS[tt] < 0 ? vv.VGIS[tt] = vv.VGIS[tt-1] + vv.Diff_VGIS[tt]*pp.f_GIS[tt] : vv.VGIS[tt-1]
           vv.T_GISstar[tt] = pp.tmaxa * (1 - vv.VGIS[tt]) ^ pp.exptstar
           vv.TD[tt] = pp.T_AT[tt] - vv.T_GISstar[tt] + 0.00001

        end

        vv.SLR_GIS[tt] = pp.slrgis * (1 - vv.VGIS[tt])

    end

end

function addGISModel(model, giscalib)

    params = CSV.read("../data/GISparams.csv", DataFrame)

    if giscalib ∉ names(params)[2:end]
        throw(ArgumentError("Unknown GIS model calibration"))
    end

    gismodel = add_comp!(model, GISModel)

    #gismodel[:meltmult] = params[params.Parameter .== "meltmult", giscalib][1]
    #gismodel[:volzero] = params[params.Parameter .== "volzero", giscalib][1]
    gismodel[:avoldot] = params[params.Parameter .== "avoldot", giscalib][1]
    gismodel[:tmaxa] = params[params.Parameter .== "tmaxa", giscalib][1]
    #gismodel[:slrgis] = params[params.Parameter .== "slrgis", giscalib][1]
    #gismodel[:expvol] = params[params.Parameter .== "expvol", giscalib][1]
    gismodel[:exptstar] = params[params.Parameter .== "exptstar", giscalib][1]

    gismodel

end
