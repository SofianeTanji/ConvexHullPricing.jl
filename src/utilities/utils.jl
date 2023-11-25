module Utilitaries
    using JSON3, DataFrames, Gurobi # Data handler
    
    const GRB_ENV = Ref{Gurobi.Env}()

    function __init__()
        GRB_ENV[] = Gurobi.Env()
        return
    end
    function ProjBox(vector, lower_bound, upper_bound)
        y = copy(vector)
        for (index, value) in enumerate(y)
            y[index] = maximum([lower_bound, minimum([value, upper_bound])])
        end
        return y
    end
    const PCU = 3000
    const PCD = -500
    const PC = abs(PCU) + abs(PCD)
    include("dataloader.jl")
    include("metrics.jl")
    include("oracles.jl")
    include("structures.jl")
end