# Implementation of the Polyak Method (see Nesterov, 2006)
using ..Utilitaries, LinearAlgebra

function PolyakMethod(instance, initial_prices, niter, ObjSol)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i=1:niter
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        PolyakStepsize = abs(fun_oracle + ObjSol)/(norm(grad_oracle))^2
        push!(iterates, iterates[i] - PolyakStepsize * grad_oracle)
    end
    push!(fun_iterates, Utilitaries.fast_oracle(instance, last(iterates))[1])
    return last(iterates), iterates, fun_iterates
end

function tPolyak(instance, opt, initial_prices, eps=1e-5)
    iterates = [initial_prices]
    fun_iterates = Array([])
    oracle_time = 0
    for i=1:500
        ((fun_oracle, grad_oracle), ortime) = @timed Utilitaries.fast_oracle(instance, iterates[i])
        oracle_time += ortime
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function

        PolyakStepsize = abs(fun_oracle + opt)/(norm(grad_oracle))^2
        push!(iterates, iterates[i] - PolyakStepsize * grad_oracle)
        if abs(fun_iterates[i] - opt)/abs(opt) < eps
            break
        end
    end
    push!(fun_iterates, Utilitaries.fast_oracle(instance, last(iterates))[1])
    return iterates, fun_iterates, oracle_time, 0
end