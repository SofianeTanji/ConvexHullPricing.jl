# Implementation of the Bundle Proximal Method (see Lemaréchal et al., 1995)
using JuMP, ..Utilitaries, LinearAlgebra, Gurobi

function BundleProximalMethod(instance, initial_prices, niter, alpha, stepsize)
    T = length(instance.Load)
    iterates = [initial_prices]
    z_iterates = [initial_prices]
    fun_iterates = []
    model_proxop = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_proxop)
    fizp1 = + Utilitaries.fast_oracle(instance, z_iterates[1])[1]
    UpperBound = Inf
    for i=1:niter
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle
        push!(fun_iterates, fun_oracle)

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end

        # Compute Z_[k + 1]
        Price = @variable(model_proxop, [1:T], lower_bound = PriceCapDown, upper_bound = PriceCapUp)
        Xi = @variable(model_proxop)
        @constraint(model_proxop, fun_oracle + dot(grad_oracle, Price - z_iterates[i])<= Xi)
        @objective(model_proxop, Min, Xi + (stepsize / 2.) * transpose(Price - iterates[i]) * (Price - iterates[i]))
        optimize!(model_proxop)
        push!(z_iterates, value.(Price))

        fzp1, gzp1 = Utilitaries.fast_oracle(instance, z_iterates[i + 1])
        fzp1, gzp1 = - fzp1, - gzp1
        fizp1 = maximum([fizp1, fzp1 + dot(gzp1, z_iterates[i + 1] - z_iterates[i])])

        if alpha * (fun_oracle - fizp1) <= fun_oracle - fzp1 # Serious Step
            push!(iterates, z_iterates[i + 1])
        else # Null step
            push!(iterates, iterates[i])
        end
    end
    return last(iterates), iterates, - fun_iterates
end

function finiteBPM(instance, initial_prices, niter, alpha, stepsize)
    T = length(instance.Load)
    iterates = [initial_prices]
    z_iterates = [initial_prices]
    fun_iterates = []
    model_proxop = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_proxop)
    fizp1 = + Utilitaries.fast_oracle(instance, z_iterates[1])[1]
    UpperBound = Inf
    for i=1:niter
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle
        push!(fun_iterates, fun_oracle)

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end

        # Compute Z_[k + 1]
        Price = @variable(model_proxop, [1:T], lower_bound = PriceCapDown, upper_bound = PriceCapUp)
        Xi = @variable(model_proxop)
        @constraint(model_proxop, fun_oracle + dot(grad_oracle, Price - z_iterates[i])<= Xi)
        @objective(model_proxop, Min, Xi + (stepsize / 2.) * transpose(Price - iterates[i]) * (Price - iterates[i]))
        optimize!(model_proxop)
        push!(z_iterates, value.(Price))
        
        fzp1, gzp1 = Utilitaries.fast_oracle(instance, z_iterates[i + 1])
        fzp1, gzp1 = - fzp1, - gzp1
        fizp1 = maximum([fizp1, fzp1 + dot(gzp1, z_iterates[i + 1] - z_iterates[i])])

        if alpha * (fun_oracle - fizp1) <= fun_oracle - fzp1 # Serious Step
            push!(iterates, z_iterates[i + 1])
        else # Null step
            push!(iterates, iterates[i])
        end
    end
    return last(iterates), iterates, - fun_iterates
end

function tBPM(instance, opt, initial_prices, alpha =.98, stepsize=1e3, eps=1e-6)
    T = length(instance.Load)
    iterates = [initial_prices]
    z_iterates = [initial_prices]
    fun_iterates = []
    model_proxop = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_proxop)
    fizp1 = + Utilitaries.fast_oracle(instance, z_iterates[1])[1]
    UpperBound = Inf
    oracle_time, solvetime = 0, 0
    for i=1:500
        ((fun_oracle, grad_oracle), ortime) = @timed Utilitaries.fast_oracle(instance, iterates[i])
        oracle_time += ortime
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle
        push!(fun_iterates, fun_oracle)

        if abs(fun_iterates[i] - opt)/abs(opt) < eps
            break
        end

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end

        # Compute Z_[k + 1]
        Price = @variable(model_proxop, [1:T], lower_bound = PriceCapDown, upper_bound = PriceCapUp)
        Xi = @variable(model_proxop)
        @constraint(model_proxop, fun_oracle + dot(grad_oracle, Price - z_iterates[i])<= Xi)
        @objective(model_proxop, Min, Xi + (stepsize / 2.) * transpose(Price - iterates[i]) * (Price - iterates[i]))
        stime = @timed optimize!(model_proxop)
        solvetime += stime.time
        push!(z_iterates, value.(Price))

        ((fzp1, gzp1), otime) = @timed Utilitaries.fast_oracle(instance, z_iterates[i + 1])
        oracle_time += otime
        fzp1, gzp1 = - fzp1, - gzp1
        fizp1 = maximum([fizp1, fzp1 + dot(gzp1, z_iterates[i + 1] - z_iterates[i])])

        if alpha * (fun_oracle - fizp1) <= fun_oracle - fzp1 # Serious Step
            push!(iterates, z_iterates[i + 1])
        else # Null step
            push!(iterates, iterates[i])
        end
    end
    return iterates, - fun_iterates, oracle_time, solvetime
end
