# Implementation of the Polyak Method (see Nesterov, 2006)
using ..Utilitaries, LinearAlgebra

function PolyakMethod(instance, initial_prices, niter, ObjSol, verbose = - 1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i=1:niter
        if verbose > 0
            @info "[Polyak: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        if (norm(grad_oracle))^2 <= 1e-5
            break
        end
        PolyakStepsize = abs(fun_oracle + ObjSol)/(norm(grad_oracle))^2
        push!(iterates, Utilitaries.ProjBox(iterates[i] - 2 * PolyakStepsize * grad_oracle), PCD, PCU)
    end
    push!(fun_iterates, Utilitaries.fast_oracle(instance, last(iterates))[1])
    return last(iterates), iterates, fun_iterates
end

function EstimatedPolyak(instance, initial_prices, niter, alpha, verbose = -1)
    iterates = [initial_prices]
    fun_iterates = Array([])

    for i=1:niter
        if verbose > 0
            @info "[Est. Polyak: Iteration $i]"
        end
        fun_oracle, grad_oracle = Utilitaries.fast_oracle(instance, iterates[i])
        push!(fun_iterates, fun_oracle)
        
        fun_oracle, grad_oracle = - fun_oracle, - grad_oracle # Maximizing a concave function <=> Minimizing a convex function
        if (norm(grad_oracle))^2 <= 1e-5
            break
        end
        PolyakStepsize = abs(fun_oracle + maximum(fun_iterates) + (alpha / sqrt(i)))/(norm(grad_oracle))^2

        push!(iterates, Utilitaries.ProjBox(iterates[i] - 2 * PolyakStepsize * grad_oracle, PCD, PCU))
    end
    x_best = iterates[argmax(fun_iterates)]
    return x_best, iterates, fun_iterates
end