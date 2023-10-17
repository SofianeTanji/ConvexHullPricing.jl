# Implementation of the Optimized Gradient Method (see Kim et al., 2016)
using ..Utilitaries, LinearAlgebra

function OptimizedGradientMethod(instance, initial_prices, niter, smoothing_parameter)
    T = length(instance.Load)
    NbGen = length(instance.ThermalGen.MinRunCapacity)
    Lips = ((1 + 1 / (NbGen^2)) * T)/smoothing_parameter
    x_iterates = [initial_prices]
    y_iterates = [initial_prices]
    Thetas = [1.]
    fun_iterates = Array([])
    for t=1:niter
        fun_oracle, grad_oracle = Utilitaries.smooth_oracle(instance, y_iterates[t], smoothing_parameter)
        push!(fun_iterates, fun_oracle)
        # fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(x_iterates, y_iterates[t] - (1 / Lips) * grad_oracle)
        if t == niter
            push!(Thetas, 0.5 * (1 + sqrt(8*Thetas[t]^2 + 1)))
        else
            push!(Thetas, 0.5 * (1 + sqrt(4*Thetas[t]^2 + 1)))
        end

        push!(y_iterates, x_iterates[t + 1] + ((Thetas[t] - 1)/Thetas[t + 1]) * (x_iterates[t + 1] - x_iterates[t]) + (Thetas[t]/Thetas[t + 1]) * (x_iterates[t + 1] - y_iterates[t]))
    end
    push!(fun_iterates, Utilitaries.fast_oracle(instance, last(x_iterates))[1])
    return last(x_iterates), x_iterates, fun_iterates
end