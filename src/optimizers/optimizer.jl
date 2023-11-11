module Optimizer

    using LinearAlgebra, Gurobi
    
    # As per current European Market Rules
    const PCU = 3000
    const PCD = -500
    const PC = abs(PCU) + abs(PCD)
    
    const GRB_ENV = Ref{Gurobi.Env}()

    function __init__()
        GRB_ENV[] = Gurobi.Env()
        return
    end

    include("dual_methods/BundleLevelMethod.jl")
    include("dual_methods/BundleProximalMethod.jl")
    include("dual_methods/DAdaptation.jl")
    include("dual_methods/DoWG.jl")
    include("dual_methods/OptimizedGradientMethod.jl")
    include("dual_methods/PolyakMethod.jl")
    include("dual_methods/SubgradientMethod.jl")

    include("primal_methods/ColumnAndRowGeneration.jl")
    include("primal_methods/ColumnGeneration.jl")
end