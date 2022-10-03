saverate = CSV.read("../data/saverate.csv", DataFrame)
saverate_world = saverate."2005-2015 average"[saverate.ISO .== "WLD"][1]

function getsaverate(iso)
    index = findfirst(saverate.ISO .== iso)

    if isnothing(index)
        0.01 * saverate_world
    else
        saverate_mine = saverate."2005-2015 average"[index]
        0.01 * (ismissing(saverate_mine) ? saverate_world : saverate_mine)
    end
end
