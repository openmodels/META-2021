@defcomp Template begin
    # Variables
    myvar = Variable(index=[time], unit="unit")

    # Parameters
    myinput = Parameter(index=[time], unit="unit")

    myparam = Parameter()

    function run_timestep(pp, vv, dd, tt)
        if is_first(tt)
            vv.myvar[tt] = pp.myinput[tt] + pp.myparam
        else
            vv.myvar[tt] = pp.myinput[tt] + vv.myvar[tt-1]
        end
    end
end

function addTemplate(model)
    template = add_comp!(model, Template)
    template[:myparam] = 100

    template
end
