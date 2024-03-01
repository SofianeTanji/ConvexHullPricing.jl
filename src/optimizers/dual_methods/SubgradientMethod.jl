# Implementation of the Subgradient Method (see Nesterov, 2006)
using ..Utilitaries, LinearAlgebra

function SubgradientMethod(instance, initial_prices, niter, alpha, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i = 1:niter
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(
            iterates,
            Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU),
        )
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end
function StochasticSubgradientMethod(
    instance,
    initial_prices,
    niter,
    alpha,
    batchsize,
    verbose = -1,
)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i = 1:niter
        if verbose > 0
            @info "[StochasticSubG: Iteration $i]"
        end
        fun_oracle, grad_oracle =
            Utilitaries.stochastic_oracle(instance, iterates[i], batchsize)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(
            iterates,
            Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU),
        )
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end
function lastSubgradientMethod(instance, initial_prices, niter, R, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for k = 1:niter
        if verbose > 0
            @info "[SubG: Iteration $k]"
        end
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[k])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        stepsize = R * (niter + 1 - k) / sqrt((niter + 1)^3)
        push!(
            iterates,
            Utilitaries.ProjBox(
                iterates[k] - stepsize * grad_oracle / norm(grad_oracle),
                PCD,
                PCU,
            ),
        )
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end

function tlastSubgradientMethod(instance, initial_prices, niter, R, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    time_vector = [0.0]
    for k = 1:niter
        if verbose > 0
            @info "[SubG: Iteration $k]"
        end
        it_time = @elapsed begin
            fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[k])
            push!(fun_iterates, fun_oracle)
            fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function

            stepsize = R * (niter + 1 - k) / sqrt((niter + 1)^3)
            push!(
                iterates,
                Utilitaries.ProjBox(
                    iterates[k] - stepsize * grad_oracle / norm(grad_oracle),
                    PCD,
                    PCU,
                ),
            )
        end
        push!(time_vector, it_time + time_vector[end])
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function tSubgradientMethod(instance, initial_prices, budget, alpha, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    time_vector = [0.0]
    i = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        it_time = @elapsed begin
            fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[i])
            push!(fun_iterates, fun_oracle)
            fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function

            push!(
                iterates,
                Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU),
            )
        end
        push!(time_vector, it_time + time_vector[end])
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function tStochasticSubgradientMethod(
    instance,
    initial_prices,
    budget,
    alpha,
    batchsize = 5,
    verbose = -1,
)
    iterates = [initial_prices]
    fun_iterates = Array([])

    time_vector = [0.0]
    i = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[StochasticSubG: Iteration $i]"
        end
        it_time = @elapsed begin
            fun_oracle, grad_oracle =
                Utilitaries.stochastic_oracle(instance, iterates[i], batchsize)
            push!(fun_iterates, fun_oracle)
            fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function

            push!(
                iterates,
                Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU),
            )
        end
        push!(time_vector, it_time + time_vector[end])
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function tRSG(instance, initial_prices, budget, K, eta1, verbose = 1)
    iterates = [initial_prices]
    fun_iterates = Array([])
    time_vector = [0.0]
    i = 1

    while time_vector[end] <= budget
        if verbose > 0
            @info "[RSG: Iteration $i]"
        end
        xb, it, fb, tb =
            tSubgradientMethod(instance, sum(iterates) / length(iterates), budget / K, eta1)
        iterates = [iterates; it]
        fun_iterates = [fun_iterates; fb]
        time_vector = [time_vector; tb[2:end] .+ time_vector[end]]
        eta1 = eta1 / 2
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function tnRSG(instance, initial_prices, budget, K, eta1, verbose = 1)
    iterates = [initial_prices]
    fun_iterates = Array([])
    time_vector = [0.0]
    i = 1

    while time_vector[end] <= budget
        if verbose > 0
            @info "[RSG: Iteration $i]"
        end
        xb, it, fb, tb = tnSubgradientMethod(
            instance,
            sum(iterates) / length(iterates),
            budget / K,
            eta1,
        )
        iterates = [iterates; it]
        fun_iterates = [fun_iterates; fb]
        time_vector = [time_vector; tb[2:end] .+ time_vector[end]]
        eta1 = eta1 / 2
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function tnSubgradientMethod(instance, initial_prices, budget, alpha, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    time_vector = [0.0]
    i = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        it_time = @elapsed begin
            fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[i])
            push!(fun_iterates, fun_oracle)
            fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function

            push!(
                iterates,
                Utilitaries.ProjBox(
                    iterates[i] - (alpha / (norm(grad_oracle) * i)) * grad_oracle,
                    PCD,
                    PCU,
                ),
            )
        end
        push!(time_vector, it_time + time_vector[end])
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function tnStochasticSubgradientMethod(
    instance,
    initial_prices,
    budget,
    alpha,
    batch,
    verbose = -1,
)
    iterates = [initial_prices]
    fun_iterates = Array([])

    time_vector = [0.0]
    i = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[Stochastic SubG: Iteration $i]"
        end
        it_time = @elapsed begin
            fun_oracle, grad_oracle =
                Utilitaries.stochastic_oracle(instance, iterates[i], batch)
            push!(fun_iterates, fun_oracle)
            fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function

            push!(
                iterates,
                Utilitaries.ProjBox(
                    iterates[i] - (alpha / (norm(grad_oracle) * i)) * grad_oracle,
                    PCD,
                    PCU,
                ),
            )
        end
        push!(time_vector, it_time + time_vector[end])
        i += 1
    end
    f_iterates = Float64[]
    for ρ in iterates
        push!(f_iterates, Utilitaries.exact_oracle(instance, ρ)[1])
    end
    x_best = iterates[argmax(f_iterates)]
    return x_best, iterates, f_iterates, time_vector
end
