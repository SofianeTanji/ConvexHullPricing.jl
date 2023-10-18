# Implementation of the DowG Method (see Khaled et al., 2023)
using ..Utilitaries, LinearAlgebra

function DoWG(instance, initial_prices, niter, initial_distance_estimate)
    list_V = [0.]
    list_R = [initial_distance_estimate]
    iterates = [initial_prices]
    fun_iterates = Array([])
    for t=1:niter
        push!(list_R, maximum([norm(iterates[t] - iterates[1]), list_R[t]]))
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[t])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        push!(list_V, list_V[t] + list_R[t + 1]^2 * norm(grad_oracle)^2)
        stepsize = (list_R[t + 1]^2)/sqrt(list_V[t + 1])
        push!(iterates, Utilitaries.ProjBox(iterates[t] - stepsize * grad_oracle, PriceCapDown, PriceCapUp))
    end
    push!(fun_iterates, Utilitaries.fast_oracle(instance, last(iterates))[1])
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end

function tDoWG(instance, opt, initial_prices, initial_distance_estimate=1e3, eps=1e-6)
    list_V = [0.]
    list_R = [initial_distance_estimate]
    iterates = [initial_prices]
    fun_iterates = Array([])
    oracle_time = 0
    for t=1:500
        push!(list_R, maximum([norm(iterates[t] - iterates[1]), list_R[t]]))
        ((fun_oracle, grad_oracle), ortime) = @timed Utilitaries.fast_oracle(instance, iterates[t])
        oracle_time += ortime
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        push!(list_V, list_V[t] + list_R[t + 1]^2 * norm(grad_oracle)^2)
        stepsize = (list_R[t + 1]^2)/sqrt(list_V[t + 1])
        push!(iterates, Utilitaries.ProjBox(iterates[t] - stepsize * grad_oracle, PriceCapDown, PriceCapUp))
        if abs(fun_iterates[t] - opt)/abs(opt) < eps
            break
        end
    end
    push!(fun_iterates, Utilitaries.fast_oracle(instance, last(iterates))[1])
    x_best = iterates[argmax(fun_iterates)]
    return iterates, fun_iterates, oracle_time, 0
end