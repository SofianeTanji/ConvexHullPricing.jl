using ..Utilitaries

function StochasticAverageGradient(instance, x0, N, α, verbose = -1)
  NbGen = length(instance.ThermalGen.MinRunCapacity)
  x_iterates = [x0]
  y = zeros(Float64, NbGen, length(x0))
  d = zeros(Float64, length(x0))
  for k=1:N
    _, gk, ik = Utilitaries.stochastic_oracle(instance, x_iterates[k], 1)
    d = d - vec(y[ik,:]) + gk
    y[ik,:] = gk
    push!(x_iterates, x_iterates[k] - (α / NbGen)*d)
  end

  # Get true values of the oracle
  f_iterates = Float64[]
  for ρ in x_iterates
        push!(f_iterates, Utilitaries.exact_oracle(instance, ρ)[1])
  end
  return last(x_iterates), x_iterates, f_iterates
end

function SAG(instance, X0, N, γ)
  T = length(instance.Load)
  G = length(instance.ThermalGen.MinRunCapacity)
  its = [X0]
  gbar = zeros(T)
  v = zeros(G, T)
  for k = 1:N
    _, g, i = Utilitaries.smooth_stochastic_oracle(instance, its[k], 1e-5, 1)
    
    gbar = gbar - vec(v[i,:])/G

    v[i,:] = g

    gbar = gbar + vec(v[i,:])/G

    push!(its, its[k] - γ * gbar)
  end
  # Get true values of the oracle
  f_iterates = Float64[]
  for ρ in its
        push!(f_iterates, Utilitaries.exact_oracle(instance, ρ)[1])
  end
  return last(its), its, f_iterates
end