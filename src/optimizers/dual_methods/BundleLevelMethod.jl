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
        fun_oracle, grad_oracle = Utilitaries.super_fast_oracle(instance, iterates[i])
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

function tOptimal(instance, initial_prices, alpha, verbose = 1)
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
    time_vector = [0.]
    idx = 1
    while UpperBound - LowerBound >= 3000
        if verbose > 0
            @info "[BLM: Iteration $idx ; Gap = $(UpperBound - LowerBound)]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[idx])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end
        @constraint(model_update_lb, fun_oracle + dot(grad_oracle, VarPrice - iterates[idx]) <= Vart)
        @objective(model_update_lb, Min, Vart)
        optimize!(model_update_lb)

        LowerBound = objective_value(model_update_lb)

        LevelSet = LowerBound + alpha * (UpperBound - LowerBound)

        ProjPrice = @variable(model_projection, [1:T], lower_bound = PCD, upper_bound = PCU)
        @constraint(model_projection, fun_oracle + dot(grad_oracle, ProjPrice - iterates[idx]) <= LevelSet)
        @objective(model_projection, Min, sum((ProjPrice[t] - iterates[idx][t])^2 for t=1:T))
        optimize!(model_projection)
        push!(iterates, value.(ProjPrice))
        end
        push!(time_vector, it_time + time_vector[end])
        idx = idx + 1
    end
    @info "UB = $(UpperBound), LB = $(LowerBound), UB-LB = $(UpperBound - LowerBound)"
    return last(iterates), iterates, fun_iterates, time_vector
end


function tBundleLevelMethod(instance, initial_prices, budget, alpha, verbose = -1)
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
    time_vector = [0.]
    idx = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[BLM: Iteration $i]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[idx])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end
        @constraint(model_update_lb, fun_oracle + dot(grad_oracle, VarPrice - iterates[idx]) <= Vart)
        @objective(model_update_lb, Min, Vart)
        optimize!(model_update_lb)

        LowerBound = objective_value(model_update_lb)

        LevelSet = LowerBound + alpha * (UpperBound - LowerBound)

        ProjPrice = @variable(model_projection, [1:T], lower_bound = PCD, upper_bound = PCU)
        @constraint(model_projection, fun_oracle + dot(grad_oracle, ProjPrice - iterates[idx]) <= LevelSet)
        @objective(model_projection, Min, sum((ProjPrice[t] - iterates[idx][t])^2 for t=1:T))
        optimize!(model_projection)
        push!(iterates, value.(ProjPrice))
        end
        push!(time_vector, it_time + time_vector[end])
        idx = idx + 1
    end
    @info "UB = $(UpperBound), LB = $(LowerBound), UB-LB = $(UpperBound - LowerBound)"
    return last(iterates), iterates, fun_iterates, time_vector
end