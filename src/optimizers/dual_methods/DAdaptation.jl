# Implementation of the D-Adaptation Method (see Defazio et al., 2023)
using ..Utilitaries, LinearAlgebra

function DAdaptation(instance, initial_prices, niter, initial_distance_estimate)
    T = length(instance.Load)
    iterates = [initial_prices]
    fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[1])
    fun_iterates = Array([])
    grad_iterates = Array([])
    S = [zeros(T)]
    D = [initial_distance_estimate]
    Gamma = [1/norm(grad_oracle)]
    for k=1:niter
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[k])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        push!(grad_iterates, grad_oracle)
        push!(S, S[k] .+ D[k] .* grad_iterates[k])
        push!(Gamma, 1/sqrt(sum(norm(grad_iterates[i])^2 for i=1:k)))
        push!(D, maximum([(Gamma[k+1]* norm(S[k+1])^2 - sum(Gamma[i]*(D[i]^2)*norm(grad_iterates[i])^2 for i=1:k) )/(2*norm(S[k+1])), D[k]]))
        push!(iterates, Utilitaries.ProjBox(iterates[1] - Gamma[k+1]*S[k+1], PriceCapDown, PriceCapUp))
    end
    coeff = 1/sum(D[j] for j=1:niter)
    final_price = coeff * sum(D[j]*iterates[j] for j=1:niter)
    push!(fun_iterates, Utilitaries.fast_oracle(instance, final_price)[1])
    return final_price, iterates, fun_iterates
end

function tDA(instance, opt, initial_prices, initial_distance_estimate=1e4, eps=1e-5)
    T = length(instance.Load)
    iterates = [initial_prices]
    fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[1])
    fun_iterates = Array([])
    grad_iterates = Array([])
    S = [zeros(T)]
    D = [initial_distance_estimate]
    Gamma = [1/norm(grad_oracle)]
    K = 500
    oracle_time = 0
    for k=1:K
        ((fun_oracle, grad_oracle), ortime) = @timed Utilitaries.fast_oracle(instance, iterates[k])
        oracle_time += ortime
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        push!(grad_iterates, grad_oracle)
        push!(S, S[k] .+ D[k] .* grad_iterates[k])
        push!(Gamma, 1/sqrt(sum(norm(grad_iterates[i])^2 for i=1:k)))
        push!(D, maximum([(Gamma[k+1]* norm(S[k+1])^2 - sum(Gamma[i]*(D[i]^2)*norm(grad_iterates[i])^2 for i=1:k) )/(2*norm(S[k+1])), D[k]]))
        push!(iterates, Utilitaries.ProjBox(iterates[1] - Gamma[k+1]*S[k+1], PriceCapDown, PriceCapUp))
        if abs(fun_iterates[k] - opt)/abs(opt) < eps
            K = k
            break
        end
    end
    coeff = 1/sum(D[j] for j=1:K)
    final_price = coeff * sum(D[j]*iterates[j] for j=1:K)
    push!(fun_iterates, Utilitaries.fast_oracle(instance, final_price)[1])
    return iterates, fun_iterates, oracle_time, 0
end