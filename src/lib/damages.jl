bhmbetas = CSV.read("../data/BHMbetas.csv", DataFrame)
amocparams = CSV.read("../data/AMOCparams.csv", DataFrame)

function getbhmbetas(iso, option, seed=nothing)
    if option == "pointestimate"
        beta1 = bhmbetas.beta1[bhmbetas.iso .== iso][1]
        beta2 = bhmbetas.beta2[bhmbetas.iso .== iso][1]
    elseif option == "low"
        beta1 = bhmbetas.beta1[bhmbetas.iso .== iso][1] - 1.96*sqrt(bhmbetas.var11[bhmbetas.iso .== iso][1])
        beta2 = bhmbetas.beta2[bhmbetas.iso .== iso][1] - 1.96*sqrt(bhmbetas.var22[bhmbetas.iso .== iso][1])
    elseif option == "high"
        beta1 = bhmbetas.beta1[bhmbetas.iso .== iso][1] + 1.96*sqrt(bhmbetas.var11[bhmbetas.iso .== iso][1])
        beta2 = bhmbetas.beta2[bhmbetas.iso .== iso][1] + 1.96*sqrt(bhmbetas.var22[bhmbetas.iso .== iso][1])
    else
        if seed != nothing
            Random.seed!(seed)
        end
        mu1, mu2 = getbhmbetas(iso, "pointestimate")
        var11 = bhmbetas.var11[bhmbetas.iso .== iso][1]
        var12 = bhmbetas.var12[bhmbetas.iso .== iso][1]
        var22 = bhmbetas.var22[bhmbetas.iso .== iso][1]
        mvn = MvNormal([mu1, mu2], [var11 var12; var12 var22])
        beta1, beta2 = rand(mvn)
    end

    return beta1, beta2
end

function gettemp1990(iso)
    return amocparams[amocparams."Country code" .== iso, "1990 temp"][1]
end

slrcoeffs = CSV.read("../data/SLRcoeffs.csv", DataFrame)

function getslrcoeff(iso, option)
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

function getslrcoeff_distribution(iso, option, qq)
    low = slrcoeffs.low[slrcoeffs.ISO .== iso][1]
    mode = slrcoeffs.mode[slrcoeffs.ISO .== iso][1]
    high = slrcoeffs.hig[slrcoeffs.ISO .== iso][1]
    quantile(TriangularDist(low, high, mode), qq)
end
