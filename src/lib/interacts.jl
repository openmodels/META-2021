function allinteractrates(callback::Function)
    params = CSV.read("../data/interactionparams.csv", DataFrame)

    outputs = Dict{String, Dict{String, Any}}()
    parammapping = Dict{String, String}("amoc" => "CMOC", "gis" => "MGIS", "wais" => "DAIS", "amaz" => "AMAZ", "nino" => "NINO")
    for forcing in ["amoc", "gis", "wais", "amaz", "nino"]
        outputs[forcing] = Dict{String, Any}()
        for affected in ["amoc", "gis", "wais", "amaz", "nino"]
            if forcing == affected
                continue
            end
            row = findfirst((params."Affected tipping point" .== parammapping[affected]) .& (params."Forcing tipping point" .== parammapping[forcing]))
            ratemu = params."Mean hazard rate factor"[row]
            ratese = params."Std. dev. of factor"[row]
            value = callback(Symbol("$(forcing)2$(affected)"), ratemu, ratese)
            outputs[forcing][affected] = value
        end
    end

    outputs
end
