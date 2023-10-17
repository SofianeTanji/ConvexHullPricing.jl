# Implementation of the Subgradient Method (see Nesterov, 2006)
using ..Utilitaries, LinearAlgebra

function SubgradientMethod(instance, initial_prices, niter, alpha)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i=1:niter
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(iterates, Utilitaries.ProjBox(iterates[i] - (alpha / sqrt(i)) * grad_oracle, PriceCapDown, PriceCapUp))
    end
    push!(fun_iterates, Utilitaries.fast_oracle(instance, last(iterates))[1])
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end
