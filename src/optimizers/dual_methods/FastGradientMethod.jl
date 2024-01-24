# Implementation of the Optimized Gradient Method (see Kim et al., 2016)
using ..Utilitaries, LinearAlgebra

function FastGradientMethod(instance, X0, niter, smoothing_parameter, verbose = -1)
    T = length(instance.Load)
    NbGen = length(instance.ThermalGen.MinRunCapacity)
    Lips = ((1 + 1 / (NbGen^2)) * T)/smoothing_parameter
    x_iterates = [X0]
    y_iterates = [X0]
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

function tFastGradientMethod(instance, X0, τ, smoothing_parameter, verbose = -1)
    T = length(instance.Load)
    NbGen = length(instance.ThermalGen.MinRunCapacity)
    Lips = 1. # ((1 + 1 / (NbGen^2)) * T)/smoothing_parameter
    x_iterates = [X0]
    y_iterates = [X0]
    fun_iterates = Array([])
    time_vector = [0.]
    t = 1
    while time_vector[end] <= τ
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

function tShiftedFGM(instance, X0, τ, α, smoothing_parameter, verbose = -1)
    x_iterates = [X0]
    y_iterates = [X0]
    shift = Utilitaries.GetShift(instance)
    fun_iterates = Array([])
    time_vector = [0.]
    t = 1
    while time_vector[end] <= τ
        if verbose > 0
            @info "[FGM: Iteration $t]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_translate_smooth_oracle(instance, y_iterates[t], smoothing_parameter, shift)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(x_iterates, Utilitaries.ProjBox(y_iterates[t] - α * grad_oracle, PCD, PCU))

        push!(y_iterates, x_iterates[t + 1] + ((t - 1)/(t + 2)) * (x_iterates[t + 1] - x_iterates[t]))
        end
        push!(time_vector, it_time + time_vector[end])
        t = t + 1
    end
    return last(x_iterates), x_iterates, fun_iterates, time_vector
end

function tGM(instance, X0, τ, α, smoothing_parameter, verbose = -1)
    iterates = [X0]
    fun_iterates = Array([])
    shift = Utilitaries.GetShift(instance)
    time_vector = [0.]
    i = 1
    while time_vector[end] <= τ
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_translate_smooth_oracle(instance, iterates[i], smoothing_parameter, shift)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(iterates, Utilitaries.ProjBox(iterates[i] - α * grad_oracle, PCD, PCU))
        end
        push!(time_vector, it_time + time_vector[end])
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function tGradientMethod(instance, X0, τ, α, smoothing_parameter, verbose = -1)
    iterates = [X0]
    fun_iterates = Array([])
    # shift = Utilitaries.GetShift(instance)
    time_vector = [0.]
    i = 1
    while time_vector[end] <= τ
        if verbose > 0
            @info "[SubG: Iteration $i]"
        end
        it_time = @elapsed begin
        fun_oracle, grad_oracle = Utilitaries.exact_smooth_oracle(instance, iterates[i], smoothing_parameter)
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        push!(iterates, Utilitaries.ProjBox(iterates[i] - α * grad_oracle, PCD, PCU))
        end
        push!(time_vector, it_time + time_vector[end])
        i += 1
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates, time_vector
end

function fista(instance, X0, Τ, μ, Lmin, L0, ρ, δ, shift)
    τ = [1. / L0]
    ρ, δ = ρ, δ
    x = [X0, X0]
    f = []
    t = 1.
    timev = [0.]
    # shift = Utilitaries.GetShift(instance)
    k = 1
    while timev[end] <= Τ
        time_it = @elapsed begin
            τk = minimum([τ[k]/δ, 1/Lmin])
            i = 0
            while true
                τn = ρ^i * τk
                tn = 0.5 * (1 + sqrt(1 + 4 * t^2 * τ[k]/τn))
                β = (t - 1)/tn
                yn = x[k + 1] + β * (x[k + 1] - x[k])
                fy, gy = Utilitaries.exact_translate_smooth_oracle(instance, yn, μ, shift)
                xn = Utilitaries.ProjBox(yn - τn * (-gy), PCD, PCU)
                fx, gx = Utilitaries.exact_smooth_oracle(instance, xn, μ)
                i += 1
                if (-fx) - (-fy) - dot(-gy, xn - yn) <= dot(xn - yn, xn - yn) / (2 * τn)
                    push!(τ, τn)
                    t = tn
                    push!(x, xn)
                    push!(f, fy)
                    break
                end
            end
        end
        push!(timev, time_it + timev[end])
        k += 1
    end
    xbest = x[2 + argmax(f)]
    return xbest, x, f, timev
end