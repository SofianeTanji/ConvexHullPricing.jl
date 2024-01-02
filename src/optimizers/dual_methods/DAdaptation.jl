# Implementation of the D-Adaptation Method (see Defazio et al., 2023)
using ..Utilitaries, LinearAlgebra

function DAdaptation(instance, initial_prices, niter, initial_distance_estimate = PC, verbose = -1)
    T = length(instance.Load)
    iterates = [initial_prices]
    fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[1])
    fun_iterates = Array([])
    grad_iterates = Array([])
    S = [zeros(T)]
    D = [initial_distance_estimate]
    Gamma = [1/norm(grad_oracle)]
    for k=1:niter
        if verbose > 0
            @info "[DAdaptation: Iteration $k]"
        end
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[k])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        push!(grad_iterates, grad_oracle)
        push!(S, S[k] .+ D[k] .* grad_iterates[k])
        push!(Gamma, 1/sqrt(sum(norm(grad_iterates[i])^2 for i=1:k)))
        push!(D, maximum([(Gamma[k+1]* norm(S[k+1])^2 - sum(Gamma[i]*(D[i]^2)*norm(grad_iterates[i])^2 for i=1:k) )/(2*norm(S[k+1])), D[k]]))
        push!(iterates, Utilitaries.ProjBox(iterates[1] - Gamma[k+1]*S[k+1], PCD, PCU))
    end
    coeff = 1/sum(D[j] for j=1:niter)
    final_price = coeff * sum(D[j]*iterates[j] for j=1:niter)
    return final_price, iterates, fun_iterates
end

function tDAdaptation(instance, initial_prices, budget, initial_distance_estimate = PC, verbose = -1)
    T = length(instance.Load)
    iterates = [initial_prices]
    fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[1])
    fun_iterates = Array([])
    grad_iterates = Array([])
    S = [zeros(T)]
    D = [initial_distance_estimate]
    Gamma = [1/norm(grad_oracle)]
    time_vector = [0.]
    k = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[DAdaptation: Iteration $k; D = $(D[end])]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[k])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        push!(grad_iterates, grad_oracle)
        push!(S, S[k] .+ D[k] .* grad_iterates[k])
        push!(Gamma, 1/sqrt(sum(norm(grad_iterates[i])^2 for i=1:k)))
        push!(D, maximum([(Gamma[k+1]* norm(S[k+1])^2 - sum(Gamma[i]*(D[i]^2)*norm(grad_iterates[i])^2 for i=1:k) )/(2*norm(S[k+1])), D[k]]))
        push!(iterates, Utilitaries.ProjBox(iterates[1] - Gamma[k+1]*S[k+1], PCD, PCU))
        end
        push!(time_vector, it_time + time_vector[end])
        k = k + 1
    end
    coeff = 1/sum(D[j] for j=1:k)
    final_price = coeff * sum(D[j]*iterates[j] for j=1:k)
    return final_price, iterates, fun_iterates, time_vector
end

function tsmoothDAdaptation(instance, initial_prices, budget, initial_distance_estimate, smoothing_parameter, verbose = -1)
    T = length(instance.Load)
    iterates = [initial_prices]
    fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[1])
    fun_iterates = Array([])
    grad_iterates = Array([])
    S = [zeros(T)]
    D = [initial_distance_estimate]
    Gamma = [1/norm(grad_oracle)]
    time_vector = [0.]
    k = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[DAdaptation: Iteration $k]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_smooth_oracle(instance, iterates[k], smoothing_parameter)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        push!(grad_iterates, grad_oracle)
        push!(S, S[k] .+ D[k] .* grad_iterates[k])
        push!(Gamma, 1/sqrt(sum(norm(grad_iterates[i])^2 for i=1:k)))
        push!(D, maximum([(Gamma[k+1]* norm(S[k+1])^2 - sum(Gamma[i]*(D[i]^2)*norm(grad_iterates[i])^2 for i=1:k) )/(2*norm(S[k+1])), D[k]]))
        push!(iterates, Utilitaries.ProjBox(iterates[1] - Gamma[k+1]*S[k+1], PCD, PCU))
        end
        push!(time_vector, it_time + time_vector[end])
        k = k + 1
    end
    coeff = 1/sum(D[j] for j=1:k)
    final_price = coeff * sum(D[j]*iterates[j] for j=1:k)
    return final_price, iterates, fun_iterates, time_vector
end