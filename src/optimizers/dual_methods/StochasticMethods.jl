using ..Utilitaries

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