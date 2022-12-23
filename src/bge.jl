using Mimi

function calculate_BGE(model::Model, emuc::Float64)
    # Present value of discounted utility streams
    sum_world_disc_utility                  = sum(model[:BGE, :world_disc_utility][1:191])
    sum_world_disc_utility_counterfactual   = sum(model[:BGE, :world_disc_utility_counterfactual][1:191])
    
    # Prefill one dimensional array of zeros for each country
    isos = dim_keys(model, :country)
    sum_disc_utility                        = zeros(length(isos), 1)
    sum_disc_utility_counterfactual         = zeros(length(isos), 1)

    for cc in 1:length(isos)
        sum_disc_utility[cc]                = sum(model[:BGE, :disc_utility][1:191, cc])
        sum_disc_utility_counterfactual[cc] = sum(model[:BGE, :disc_utility_counterfactual][1:191, cc])
    end

    # Apply BGE method (eq. 5 in Anthoff and Tol 2009 ERE)
    world_bge = (sum_world_disc_utility_counterfactual/sum_world_disc_utility)^(1/(1-emuc))-1 
    
    #THIS IS PROBABLY CORRECT; BUT NEED TO FIGURE OUT HOW TO HAVE THE FUNCTION PRODUCE MORE THAN JUST ONE SCALAR. THAT SEEMS TO BE THE CURRENT PROBLEM
    #bge                                     = zeros(length(isos), 1)
    #for cc in 1:length(isos)
    #   bge[cc] = (sum_disc_utility_counterfactual[cc]/sum_disc_utility[cc])^(1/(1-emuc))-1
    #end
end

include("../src/MimiMETA.jl")
model = base_model(; rcp="CP-GMP", tdamage="pointestimate", slrdamage="mode")
run(model)
calculate_BGE(model, 1.5)
