# Implementation of the Bundle Proximal Method (see Lemaréchal et al., 1995)
using JuMP, ..Utilitaries, LinearAlgebra, Gurobi

mPCD = -500.
mPCU = 3000.

function BundleProximalMethod(instance, initial_prices, niter, α, stepsize, verbose = -1)
    T = length(instance.Load)
    iterates = [initial_prices]
    z_iterates = [initial_prices]
    fun_iterates = []
    model_proxop = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_proxop)
    set_optimizer_attributes(model_proxop, "MIPGap" => 0, "MIPGapAbs" => 1e-8)

    fizp1 = + Utilitaries.exact_oracle(instance, z_iterates[1])[1]
    UpperBound = Inf
    for i=1:niter
        if verbose > 0
            @info "[BPM: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[i])
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle
        push!(fun_iterates, fun_oracle)

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end

        # Compute Z_[k + 1]
        Price = @variable(model_proxop, [1:T], lower_bound = mPCD, upper_bound = mPCU)
        Xi = @variable(model_proxop)
        @constraint(model_proxop, fun_oracle + dot(grad_oracle, Price - z_iterates[i])<= Xi)
        @objective(model_proxop, Min, Xi + (stepsize / 2.) * transpose(Price - iterates[i]) * (Price - iterates[i]))
        optimize!(model_proxop)
        push!(z_iterates, value.(Price))

        fzp1, gzp1 = Utilitaries.exact_oracle(instance, z_iterates[i + 1])
        fzp1, gzp1 = - fzp1, - gzp1
        fizp1 = maximum([fizp1, fzp1 + dot(gzp1, z_iterates[i + 1] - z_iterates[i])])

        if α * (fun_oracle - fizp1) <= fun_oracle - fzp1 # Serious Step
            push!(iterates, z_iterates[i + 1])
        else # Null step
            push!(iterates, iterates[i])
        end
    end
    return last(iterates), iterates, - fun_iterates
end

function tBundleProximalMethod(instance, initial_prices, τ, α, stepsize, verbose = -1)
    T = length(instance.Load)
    iterates = [initial_prices]
    z_iterates = [initial_prices]
    fun_iterates = []
    model_proxop = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_proxop)
    set_optimizer_attributes(model_proxop, "MIPGap" => 0, "MIPGapAbs" => 1e-8)

    fizp1 = + Utilitaries.exact_oracle(instance, z_iterates[1])[1]
    UpperBound = Inf
    time_vector = [0.]
    idx = 1
    while time_vector[end] <= τ
        if verbose > 0
            @info "[BPM: Iteration $(idx)]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[idx])
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle
        push!(fun_iterates, fun_oracle)

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end

        # Compute Z_[k + 1]
        Price = @variable(model_proxop, [1:T], lower_bound = mPCD, upper_bound = mPCU)
        Xi = @variable(model_proxop)
        @constraint(model_proxop, fun_oracle + dot(grad_oracle, Price - z_iterates[idx])<= Xi)
        @objective(model_proxop, Min, Xi + (stepsize / 2.) * transpose(Price - iterates[idx]) * (Price - iterates[idx]))
        optimize!(model_proxop)
        push!(z_iterates, value.(Price))

        fzp1, gzp1 = Utilitaries.exact_oracle(instance, z_iterates[idx + 1])
        fzp1, gzp1 = - fzp1, - gzp1
        fizp1 = maximum([fizp1, fzp1 + dot(gzp1, z_iterates[idx + 1] - z_iterates[idx])])

        if α * (fun_oracle - fizp1) <= fun_oracle - fzp1 # Serious Step
            push!(iterates, z_iterates[idx + 1])
        else # Null step
            push!(iterates, iterates[idx])
        end
        end
        push!(time_vector, it_time + time_vector[end])
        idx = idx + 1
    end
    return last(iterates), iterates, -fun_iterates, time_vector
end

function BundleProximalLevelMethod(instance, initial_prices, niter, α, verbose = -1)

    T = length(instance.Load)
    iterates = [initial_prices]
    fun_iterates = Array([])

    UpperBound, LowerBound = Inf, - Inf

    model_update_lb = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_update_lb)
    set_optimizer_attributes(model_update_lb, "MIPGap" => 0, "MIPGapAbs" => 1e-8)


    VarPrice = @variable(model_update_lb, [1:T], lower_bound = mPCD, upper_bound = mPCU)
    Vart = @variable(model_update_lb)
    
    model_projection = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_projection)
    set_optimizer_attributes(model_projection, "MIPGap" => 0, "MIPGapAbs" => 1e-8)

    newGap = Inf
    newLevel = Inf
    BestPoint = initial_prices
    prevFunVal = Inf

    for i=1:niter
        if verbose > 0
            @info "[BPLM: Iteration $i]"
        end

        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        if fun_oracle < prevFunVal
            BestPoint = iterates[i]
        end

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end

        @constraint(model_update_lb, fun_oracle + dot(grad_oracle, VarPrice - iterates[i]) <= Vart)
        @objective(model_update_lb, Min, Vart)
        optimize!(model_update_lb)

        LowerBound = objective_value(model_update_lb)

        LevelSet = LowerBound + α * (UpperBound - LowerBound)

        if UpperBound - LowerBound >= (1 - α) * newGap
            newLevel = minimum([LevelSet, newLevel])
        else
            newLevel = LevelSet
            newGap = UpperBound - LowerBound
        end

        ProjPrice = @variable(model_projection, [1:T], lower_bound = mPCD, upper_bound = mPCU)
        @constraint(model_projection, fun_oracle + dot(grad_oracle, ProjPrice - iterates[i]) <= newLevel)
        @objective(model_projection, Min, sum((ProjPrice[t] - BestPoint[t])^2 for t=1:T))
        optimize!(model_projection)
        push!(iterates, value.(ProjPrice))
        prevFunVal = fun_oracle
    end
    return last(iterates), iterates, fun_iterates
end

function tBundleProximalLevelMethod(instance, initial_prices, τ, α, verbose = -1)

    T = length(instance.Load)
    iterates = [initial_prices]
    fun_iterates = Array([])

    UpperBound, LowerBound = Inf, - Inf
    UpperBounds = Float64[]
    LowerBounds = Float64[]

    model_update_lb = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_update_lb)
    set_optimizer_attributes(model_update_lb, "MIPGap" => 0, "MIPGapAbs" => 1e-8)


    VarPrice = @variable(model_update_lb, [1:T], lower_bound = mPCD, upper_bound = mPCU)
    Vart = @variable(model_update_lb)
    
    model_projection = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_projection)
    set_optimizer_attributes(model_projection, "MIPGap" => 0, "MIPGapAbs" => 1e-8)

    newGap = Inf
    newLevel = Inf
    BestPoint = initial_prices
    prevFunVal = Inf

    time_vector = [0.]
    i = 1
    while time_vector[end] <= τ
        if verbose > 0
            @info "[BPLM: Iteration $(i) ; UB = $(UpperBound), LB = $(LowerBound), UB-LB = $(UpperBound - LowerBound)]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        if fun_oracle < prevFunVal
            BestPoint = iterates[i]
        end

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end
        push!(UpperBounds, UpperBound)

        @constraint(model_update_lb, fun_oracle + dot(grad_oracle, VarPrice - iterates[i]) <= Vart)
        @objective(model_update_lb, Min, Vart)
        optimize!(model_update_lb)

        LowerBound = objective_value(model_update_lb)
        push!(LowerBounds, LowerBound)

        LevelSet = LowerBound + α * (UpperBound - LowerBound)

        if UpperBound - LowerBound >= (1- α) * newGap
            newLevel = minimum([LevelSet, newLevel])
        else
            newLevel = LevelSet
            newGap = UpperBound - LowerBound
        end

        ProjPrice = @variable(model_projection, [1:T], lower_bound = mPCD, upper_bound = mPCU)
        @constraint(model_projection, fun_oracle + dot(grad_oracle, ProjPrice - iterates[i]) <= newLevel)
        @objective(model_projection, Min, sum((ProjPrice[t] - BestPoint[t])^2 for t=1:T))
        optimize!(model_projection)
        push!(iterates, value.(ProjPrice))
        prevFunVal = fun_oracle
        end
        push!(time_vector, it_time + time_vector[end])
        i = i + 1
    end
    @info "UB = $(UpperBound), LB = $(LowerBound), UB-LB = $(UpperBound - LowerBound)"
    return last(iterates), iterates, fun_iterates, time_vector, UpperBounds, LowerBounds
end

function tBPLM(instance, initial_prices, τ, α, verbose = -1)

    T = length(instance.Load)
    iterates = [initial_prices]
    fun_iterates = Array([])

    UpperBound, LowerBound = Inf, - Inf

    model_update_lb = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_update_lb)
    set_optimizer_attributes(model_update_lb, "MIPGap" => 0, "MIPGapAbs" => 1e-8)


    VarPrice = @variable(model_update_lb, [1:T], lower_bound = mPCD, upper_bound = mPCU)
    Vart = @variable(model_update_lb)
    
    model_projection = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model_projection)
    set_optimizer_attributes(model_projection, "MIPGap" => 0, "MIPGapAbs" => 1e-8)

    newGap = Inf
    newLevel = Inf
    BestPoint = initial_prices
    prevFunVal = Inf

    time_vector = [0.]
    i = 1
    while time_vector[end] <= τ
        if verbose > 0
            @info "[BPLM: Iteration $(i) ; UB = $(UpperBound), LB = $(LowerBound), UB-LB = $(UpperBound - LowerBound)]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        if fun_oracle < prevFunVal
            BestPoint = iterates[i]
        end

        if UpperBound > fun_oracle
            UpperBound = fun_oracle
        end

        @constraint(model_update_lb, fun_oracle + dot(grad_oracle, VarPrice - iterates[i]) <= Vart)
        @objective(model_update_lb, Min, Vart)
        optimize!(model_update_lb)

        LowerBound = objective_value(model_update_lb)

        LevelSet = LowerBound + α * (UpperBound - LowerBound)

        if UpperBound - LowerBound >= (1- α) * newGap
            newLevel = minimum([LevelSet, newLevel])
        else
            newLevel = LevelSet
            newGap = UpperBound - LowerBound
        end

        ProjPrice = @variable(model_projection, [1:T], lower_bound = mPCD, upper_bound = mPCU)
        @constraint(model_projection, fun_oracle + dot(grad_oracle, ProjPrice - iterates[i]) <= newLevel)
        @objective(model_projection, Min, sum((ProjPrice[t] - iterates[i][t])^2 for t=1:T))
        optimize!(model_projection)
        push!(iterates, value.(ProjPrice))
        prevFunVal = fun_oracle
        end
        push!(time_vector, it_time + time_vector[end])
        i = i + 1
    end
    @info "UB = $(UpperBound), LB = $(LowerBound), UB-LB = $(UpperBound - LowerBound)"
    return last(iterates), iterates, fun_iterates, time_vector
end