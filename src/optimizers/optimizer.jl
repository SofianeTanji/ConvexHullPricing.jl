module Optimizer
using LinearAlgebra, JuMP, Gurobi

const GRB_ENV = Ref{Gurobi.Env}()

function __init__()
    GRB_ENV[] = Gurobi.Env()
    return
end

const PCU = 300.0
const PCD = -300.0
const PC = abs(PCU) + abs(PCD)
include("dual_methods/BundleLevelMethod.jl")
include("dual_methods/BundleProximalMethod.jl")
include("dual_methods/DAdaptation.jl")
include("dual_methods/DoWG.jl")
include("dual_methods/FastGradientMethod.jl")
include("dual_methods/PolyakMethod.jl")
include("dual_methods/SubgradientMethod.jl")

include("primal_methods/ColumnAndRowGeneration.jl")
include("primal_methods/ColumnGeneration.jl")
end
