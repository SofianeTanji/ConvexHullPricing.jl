# Implementation of the Bundle Level Method (see Lemaréchal et al., 1995)
using JuMP, ..Utilitaries, LinearAlgebra, Gurobi

function BundleLevelMethod(instance, initial_prices, niter, alpha, verbose = -1)
    T = length(instance.Load)
    iterates = [initial_prices]
    fun_iterates = Array([])

    UpperBound, LowerBound = Inf, - Inf

    model_update_lb = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_update_lb)

    VarPrice = @variable(model_update_lb, [1:T], lower_bound = PCD, upper_bound = PCU)
    Vart = @variable(model_update_lb)
    
    model_projection = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_projection)

    for i=1:niter
        if verbose > 0
            @info "[BLM: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end
        @constraint(model_update_lb, fun_oracle + dot(grad_oracle, VarPrice - iterates[i]) <= Vart)
        @objective(model_update_lb, Min, Vart)
        optimize!(model_update_lb)

        LowerBound = objective_value(model_update_lb)

        LevelSet = LowerBound + alpha * (UpperBound - LowerBound)

        ProjPrice = @variable(model_projection, [1:T], lower_bound = PCD, upper_bound = PCU)
        @constraint(model_projection, fun_oracle + dot(grad_oracle, ProjPrice - iterates[i]) <= LevelSet)
        @objective(model_projection, Min, sum((ProjPrice[t] - iterates[i][t])^2 for t=1:T))
        optimize!(model_projection)
        push!(iterates, value.(ProjPrice))
    end
    return last(iterates), iterates, fun_iterates
end

function GetOptimal(instance, initial_prices, itermax = 100, RelativeGap = 1e-7, alpha = .98, verbose = -1)
    T = length(instance.Load)
    iterates = [initial_prices]
    fun_iterates = Array([])

    UpperBound, LowerBound, RelGap = Inf, - Inf, Inf

    model_update_lb = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_update_lb)

    VarPrice = @variable(model_update_lb, [1:T], lower_bound = PCD, upper_bound = PCU)
    Vart = @variable(model_update_lb)
    
    model_projection = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_projection)

    for i=1:itermax
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end
        @constraint(model_update_lb, fun_oracle + dot(grad_oracle, VarPrice - iterates[i]) <= Vart)
        @objective(model_update_lb, Min, Vart)
        optimize!(model_update_lb)

        LowerBound = objective_value(model_update_lb)

        RelGap = abs(UpperBound - LowerBound)/maximum([abs(UpperBound), abs(LowerBound)])
        if verbose > 0
            @show "Current iteration is $i. Relgap is $RelGap"
        end
        if RelGap <= RelativeGap
            break
        end

        LevelSet = LowerBound + alpha * (UpperBound - LowerBound)

        ProjPrice = @variable(model_projection, [1:T], lower_bound = PCD, upper_bound = PCU)
        @constraint(model_projection, fun_oracle + dot(grad_oracle, ProjPrice - iterates[i]) <= LevelSet)
        @objective(model_projection, Min, sum((ProjPrice[t] - iterates[i][t])^2 for t=1:T))
        optimize!(model_projection)
        push!(iterates, value.(ProjPrice))
    end
    push!(fun_iterates, Utilitaries.fast_oracle(instance, last(iterates))[1])
    return last(iterates), iterates, fun_iterates, RelGap
end