# Implementation of the Subgradient Method (see Nesterov, 2006)
using ..Utilitaries, LinearAlgebra

function SubgradientMethod(instance, initial_prices, niter, alpha, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i=1:niter
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.super_fast_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(iterates, Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU))
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end
function StochasticSubgradientMethod(instance, initial_prices, niter, alpha, batchsize, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i=1:niter
        if verbose > 0
            @info "[StochasticSubG: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.stochastic_oracle(instance, iterates[i], batchsize)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(iterates, Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU))
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end
function lastSubgradientMethod(instance, initial_prices, niter, R, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for k=1:niter
        if verbose > 0
            @info "[SubG: Iteration $k]"
        end
        fun_oracle, grad_oracle = Utilitaries.super_fast_oracle(instance, iterates[k])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        
        stepsize = R * (niter + 1 - k)/sqrt((niter+1)^3)
        push!(iterates, Utilitaries.ProjBox(iterates[k] - stepsize * grad_oracle / norm(grad_oracle), PCD, PCU))
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end


function tSubgradientMethod(instance, initial_prices, budget, alpha, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    time_vector = [0.]
    i = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(iterates, Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU))
        end
        push!(time_vector, it_time + time_vector[end])
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end
function tStochasticSubgradientMethod(instance, initial_prices, budget, alpha, batchsize, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    time_vector = [0.]
    i = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[StochasticSubG: Iteration $i]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.stochastic_oracle(instance, iterates[i], batchsize)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(iterates, Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU))
        end
        push!(time_vector, it_time + time_vector[end])
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end