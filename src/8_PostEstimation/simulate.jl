###############################################################################################
# Simulate set of models for shocks passed to function
###############################################################################################
function simulate(sr, lr, init_conds, shocks, meas_err, select_variables, simnames, timeline)

    T = size(shocks[1], 2)
    n_shocks = length(shocks)
    n_select_var = length(select_variables)
    SIMs = Array{Array{Float64}}(undef, n_shocks)
    for j = 1:n_shocks
        selector = []
        obs_selector = []
        obs_sel = []
        iter = 1
        for i in select_variables
            try
                idx = getfield(sr.indexes_r, i)
                append!(selector, idx)
                if i in e_set.observed_vars_input
                    append!(obs_selector, findfirst(==(i),e_set.observed_vars_input))
                    append!(obs_sel, iter)
                end
            catch
                append!(selector, sr.n_par.ntotal_r + 1)
                println(selector)
            end
            iter += 1
        end

        SIMs_aux = zeros(length(selector), T)
        x = init_conds[:, 1]
        MX = [I; lr.State2Control; zeros(1, sr.n_par.nstates_r)]
        for t = 1:T
            SIMs_aux[:, t] = MX[selector, :] * x
            x[:] = lr.LOMstate * x
            x[:] += shocks[j][:, t] # shock in "t" moves observables in "t+1" 
            SIMs_aux[obs_sel,t] += meas_err[obs_selector,t]
        end
        SIMs[j] = SIMs_aux
        data = vcat(
            [
                DataFrame(
                    Time = timeline[t],
                    Variable = select_variables[v],
                    Value = SIMs_aux[v, t],
                ) for t = 1:T, v = 1:n_select_var
            ]...,
            )
        CSV.write(string("8_PostEstimation/Tables/Sim_",simnames[j],".csv"),
        unstack(data, :Variable, :Value)
        )
    end
    return SIMs

end
