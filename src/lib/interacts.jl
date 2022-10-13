function allinteractcalls(callback::Function)
    outputs = Dict{String, Dict{String, Any}}()
    parammapping = Dict{String, String}("amoc" => "CMOC", "gis" => "MGIS", "wais" => "DAIS", "amaz" => "AMAZ", "nino" => "NINO")
    for forcing in ["amoc", "gis", "wais", "amaz", "nino"]
        outputs[forcing] = Dict{String, Any}()
        for affected in ["amoc", "gis", "wais", "amaz", "nino"]
            if forcing == affected
                continue
            end
            outputs[forcing][affected] = callback(Symbol("$(forcing)2$(affected)"), parammapping[affected], parammapping[forcing])
        end
    end

    outputs
end    

function allinteractrates(callback::Function)
    params = CSV.read("../data/interactionparams.csv", DataFrame)

    function innercallback(symb, paramaffected, paramforcing)
        row = findfirst((params."Affected tipping point" .== paramaffected) .& (params."Forcing tipping point" .== paramforcing))
        ratemu = params."Mean hazard rate factor"[row]
        ratese = params."Std. dev. of factor"[row]
        value = callback(symb, ratemu, ratese)
    end
    
    allinteractcalls(innercallback)
end
