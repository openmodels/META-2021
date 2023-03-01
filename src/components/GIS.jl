@defcomp GISModel begin

    # Variables

    VGIS = Variable(index=[time], unit="percent") # volume of the GIS relative to initial volume
    Diff_VGIS = Variable(index=[time], unit="percent") # difference in VGIS between periods
    T_GISstar = Variable(index=[time], unit="degC") # temperature increase relative to initial temperature associated with a particular degree of GIS melting in equilibrium
    TD = Variable(index=[time], unit="degC") # difference between T_AT and T_GISstar
    SLR_GIS = Variable(index=[time], unit="m") # sea-level rise from GIS melting

    avoldot = Variable() # Equilibrium volume temperature equation

    # Parameters

    T_AT = Parameter(index=[time], unit="degC") # Atmospheric temperature relative to pre-industrial, output of temperature model

    meltmult = Parameter(default=1) # Multiplier times melt rate
    volzero = Parameter(default=0.999, unit="percent") # Initial volume of GIS (%)
    avoldot0 = Parameter() # Nordhaus quantifies s.e. of this parameter as 2.44E-05 (/5/100). This comes out of regressing the equation determining T*(t) on data from Robinson et al. (2012). Should be used with tmaxa=3.4 only. Nordhaus also carries out sensitivity analysis with double melt rate but says this is beyond any reported estimates in AR5. Note low and high scenarios defined as 2.5% and 97.5% of distribution.
    tmaxa = Parameter(unit="degC") # Global temperature at which minimum volume
    slrgis = Parameter(default=7, unit="metres") # Sea level rise from GIS melt (metres)
    expvol = Parameter(default=0.2) # Exponent on voldot
    exptstar = Parameter() # Exponent on Tstar

    f_GIS = Parameter(index=[time])

    function init(pp, vv, dd)
        vv.avoldot = pp.avoldot0 * pp.meltmult
    end

    function run_timestep(pp, vv, dd, tt)

        if is_first(tt)

           vv.VGIS[tt] = pp.volzero
           vv.T_GISstar[tt] = pp.tmaxa * (1 - vv.VGIS[tt]) ^ pp.exptstar
           vv.TD[tt] = pp.T_AT[tt] - vv.T_GISstar[tt] + 0.00001

        else

           vv.Diff_VGIS[tt] = vv.avoldot * sign(vv.TD[tt-1]) * vv.TD[tt-1] ^ 2 * vv.VGIS[tt-1] ^ pp.expvol
           vv.VGIS[tt] = (vv.Diff_VGIS[tt] < 0 ? vv.VGIS[tt-1] + vv.Diff_VGIS[tt]*pp.f_GIS[tt] : vv.VGIS[tt-1])
           vv.T_GISstar[tt] = pp.tmaxa * (1 - vv.VGIS[tt]) ^ pp.exptstar
           vv.TD[tt] = pp.T_AT[tt] - vv.T_GISstar[tt] + 0.00001

        end

        vv.SLR_GIS[tt] = pp.slrgis * (1 - vv.VGIS[tt])

    end

end

function addGISModel(model, giscalib; before=nothing, after=nothing)

    params = CSV.read("../data/GISparams.csv", DataFrame)

    if giscalib ∉ names(params)[2:end]
        throw(ArgumentError("Unknown GIS model calibration"))
    end

    gismodel = add_comp!(model, GISModel, before=before, after=after)

    #gismodel[:meltmult] = params[params.Parameter .== "meltmult", giscalib][1]
    #gismodel[:volzero] = params[params.Parameter .== "volzero", giscalib][1]
    gismodel[:avoldot0] = params[params.Parameter .== "avoldot0", giscalib][1]
    gismodel[:tmaxa] = params[params.Parameter .== "tmaxa", giscalib][1]
    #gismodel[:slrgis] = params[params.Parameter .== "slrgis", giscalib][1]
    #gismodel[:expvol] = params[params.Parameter .== "expvol", giscalib][1]
    gismodel[:exptstar] = params[params.Parameter .== "exptstar", giscalib][1]

    gismodel[:f_GIS] = ones(dim_count(model, :time))

    gismodel

end
