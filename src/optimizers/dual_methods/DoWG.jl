# Implementation of the DowG Method (see Khaled et al., 2023)
using ..Utilitaries, LinearAlgebra

function DowG(instance, initial_prices, niter, initial_distance_estimate = PC, verbose = -1)
    list_V = [0.0]
    list_R = [initial_distance_estimate]
    iterates = [initial_prices]
    fun_iterates = Array([])
    for t = 1:niter
        if verbose > 0
            @info "[DowG: Iteration $t]"
        end
        push!(list_R, maximum([norm(iterates[t] - iterates[1]), list_R[t]]))
        fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[t])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        push!(list_V, list_V[t] + list_R[t+1]^2 * norm(grad_oracle)^2)
        stepsize = (list_R[t+1]^2) / sqrt(list_V[t+1])
        push!(iterates, Utilitaries.ProjBox(iterates[t] - stepsize * grad_oracle, PCD, PCU))
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end


function tDowG(
    instance,
    initial_prices,
    budget,
    initial_distance_estimate = PC,
    verbose = -1,
)
    list_V = [0.0]
    list_R = [initial_distance_estimate]
    iterates = [initial_prices]
    fun_iterates = Array([])
    time_vector = [0.0]
    t = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[DowG: Iteration $t]"
        end
        it_time = @elapsed begin
            push!(list_R, maximum([norm(iterates[t] - iterates[1]), list_R[t]]))
            fun_oracle, grad_oracle = Utilitaries.exact_oracle(instance, iterates[t])
            push!(fun_iterates, fun_oracle)
            fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function
            push!(list_V, list_V[t] + list_R[t+1]^2 * norm(grad_oracle)^2)
            stepsize = (list_R[t+1]^2) / sqrt(list_V[t+1])
            push!(
                iterates,
                Utilitaries.ProjBox(iterates[t] - stepsize * grad_oracle, PCD, PCU),
            )
        end
        push!(time_vector, it_time + time_vector[end])
        t = t + 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function tsmoothDowG(
    instance,
    initial_prices,
    budget,
    initial_distance_estimate,
    smoothing_parameter,
    verbose = -1,
)
    list_V = [0.0]
    list_R = [initial_distance_estimate]
    iterates = [initial_prices]
    fun_iterates = Array([])
    time_vector = [0.0]
    t = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[DowG: Iteration $t]"
        end
        it_time = @elapsed begin
            push!(list_R, maximum([norm(iterates[t] - iterates[1]), list_R[t]]))
            fun_oracle, grad_oracle = Utilitaries.exact_smooth_oracle(
                instance,
                iterates[t],
                smoothing_parameter,
            )
            push!(fun_iterates, fun_oracle)
            fun_oracle, grad_oracle = -fun_oracle, -grad_oracle # Maximizing a concave function <=> Minimizing a convex function
            push!(list_V, list_V[t] + list_R[t+1]^2 * norm(grad_oracle)^2)
            stepsize = (list_R[t+1]^2) / sqrt(list_V[t+1])
            push!(
                iterates,
                Utilitaries.ProjBox(iterates[t] - stepsize * grad_oracle, PCD, PCU),
            )
        end
        push!(time_vector, it_time + time_vector[end])
        t = t + 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end
