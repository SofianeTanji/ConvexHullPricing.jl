# Implementation of the Optimized Gradient Method (see Kim et al., 2016)
using ..Utilitaries, LinearAlgebra

function FastGradientMethod(instance, initial_prices, niter, smoothing_parameter, verbose = -1)
    T = length(instance.Load)
    NbGen = length(instance.ThermalGen.MinRunCapacity)
    Lips = ((1 + 1 / (NbGen^2)) * T)/smoothing_parameter
    x_iterates = [initial_prices]
    y_iterates = [initial_prices]
    fun_iterates = Array([])
    for t=1:niter
        if verbose > 0
            @info "[FGM: Iteration $t]"
        end    
        fun_oracle, grad_oracle = Utilitaries.exact_smooth_oracle(instance, y_iterates[t], smoothing_parameter)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(x_iterates, y_iterates[t] - (1 / Lips) * grad_oracle)

        push!(y_iterates, x_iterates[t + 1] + ((t - 1)/(t + 2)) * (x_iterates[t + 1] - x_iterates[t]))
    end
    return last(x_iterates), x_iterates, fun_iterates
end

function tFastGradientMethod(instance, initial_prices, budget, smoothing_parameter, verbose = -1)
    T = length(instance.Load)
    NbGen = length(instance.ThermalGen.MinRunCapacity)
    Lips = ((1 + 1 / (NbGen^2)) * T)/smoothing_parameter
    x_iterates = [initial_prices]
    y_iterates = [initial_prices]
    fun_iterates = Array([])
    time_vector = [0.]
    t = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[FGM: Iteration $t]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_smooth_oracle(instance, y_iterates[t], smoothing_parameter)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(x_iterates, y_iterates[t] - (1 / Lips) * grad_oracle)

        push!(y_iterates, x_iterates[t + 1] + ((t - 1)/(t + 2)) * (x_iterates[t + 1] - x_iterates[t]))
        end
        push!(time_vector, it_time + time_vector[end])
        t = t + 1
    end
    return last(x_iterates), x_iterates, fun_iterates, time_vector
end