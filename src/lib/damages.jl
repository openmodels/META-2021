import Random

bhmbetas = CSV.read("../data/BHMbetas.csv", DataFrame)
amocparams = CSV.read("../data/AMOCparams.csv", DataFrame)

function getbhmbetas(tdamage::String)
    if tdamage == "pointestimate"
        return 0.0127183, -0.0004871
    elseif tdamage == "low"
        beta2 = -0.0004871 + 1.96 * 0.0001183216 # beta2 + 1.96 SE_2
        beta1 = -2 * beta2 * 13.05512 # keep same peak
        return beta1, beta2
    elseif tdamage == "high"
        beta2 = -0.0004871 - 1.96 * 0.0001183216 # beta2 + 1.96 SE_2
        beta1 = -2 * beta2 * 13.05512 # keep same peak
        return beta1, beta2
    else # distribution
        mvn = Distributions.MvNormal([0.0127183, -0.0004871], [1.43e-5 -3.76e-7; -3.76e-7 1.4e-8])
        rand(mvn)
    end
end

function gettemp1990(iso)
    return amocparams[amocparams."Country code" .== iso, "1990 temp"][1]
end

slrcoeffs = CSV.read("../data/SLRcoeffs.csv", DataFrame)

function getslrcoeff(iso, option::String)
    if option == "mode"
        return slrcoeffs.mode[slrcoeffs.ISO .== iso][1]
    elseif option == "low"
        return slrcoeffs.low[slrcoeffs.ISO .== iso][1]
    elseif option == "high"
        return slrcoeffs.hig[slrcoeffs.ISO .== iso][1]
    else
        throw(ArgumentError("Distribution implemented in getslrcoeff_distribution"))
    end
end

function getslrcoeff_distribution(iso, option::String, qq)
    low = slrcoeffs.low[slrcoeffs.ISO .== iso][1]
    mode = slrcoeffs.mode[slrcoeffs.ISO .== iso][1]
    high = slrcoeffs.hig[slrcoeffs.ISO .== iso][1]
    quantile(TriangularDist(low, high, mode), qq)
end
