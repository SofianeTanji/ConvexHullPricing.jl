# Implementation of the Subgradient Method (see Nesterov, 2006)
using ..Utilitaries, LinearAlgebra

function SubgradientMethod(instance, initial_prices, niter, alpha, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i=1:niter
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(iterates, Utilitaries.ProjBox(iterates[i] - (alpha / (i)) * grad_oracle, PCD, PCU))
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end
function StochasticSubgradientMethod(instance, initial_prices, niter, alpha, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i=1:niter
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.stochastic_oracle(instance, iterates[i], 10)
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
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[k])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        
        stepsize = R * (niter + 1 - k)/sqrt((niter+1)^3)
        push!(iterates, Utilitaries.ProjBox(iterates[k] - stepsize * grad_oracle / norm(grad_oracle), PCD, PCU))
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end