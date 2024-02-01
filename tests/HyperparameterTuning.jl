# Load packages
using Revise
using ConvexHullPricing
using DataFrames
using Plots
using JLD2
using ProgressBars, JuMP, Gurobi

const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer

# Load instances
BEinstances = []
for file in readdir("data//belgian"; join=true)
  push!(BEinstances, UT.load_data(file))
end
CAinstances = []
for file in readdir("data//ca"; join=true)
  push!(CAinstances, UT.load_data(file))
end
instanceBE = BEinstances[8]
instanceCA = CAinstances[1]
X0BE = UT.LP_Relaxation(instanceBE)
X0CA = UT.LP_Relaxation(instanceCA)
Τ = 3 * 60.

# Method 1 : D-Adaptation
function HPO_Coarse_DA()
  CoarseParameterRange = [10.0^k for k = -4:4]
  CoarseResultsBE = []
  for α in CoarseParameterRange
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDAdaptation(instanceBE, X0BE, Τ, α)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseDA-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  for α in CoarseParameterRange
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDAdaptation(instanceCA, X0CA, Τ, α)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseDA-CA.jld2", CoarseResultsCA)
end
# HPO_Coarse_DA()

function HPO_fine_DA()
  # Saved value : α = 10
  FineParameterRangeBE = [5., 10., 15., 20., 25.]
  FineResultsBE = []
  for α in FineParameterRangeBE
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDAdaptation(instanceBE, X0BE, Τ, α)
    push!(FineResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//FineDA-BE.jld2", FineResultsBE)
  # Saved value : α = 0.1
  FineParameterRangeCA = [0.05, 0.1, 0.15]
  FineResultsCA = []
  for α in FineParameterRangeCA
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDAdaptation(instanceCA, X0CA, Τ, α)
    push!(FineResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//FineDA-CA.jld2", FineResultsCA)
end
# HPO_fine_DA()

# Method 1 : DOWG
function HPO_Coarse_DOWG()
  CoarseParameterRange = [10.0^k for k = -4:4]
  CoarseResultsBE = []
  for α in CoarseParameterRange
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDowG(instanceBE, X0BE, Τ, α)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseDowG-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  for α in CoarseParameterRange
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDowG(instanceCA, X0CA, Τ, α)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseDowG-CA.jld2", CoarseResultsCA)
end
# HPO_Coarse_DA()

function HPO_Coarse_SUBG()
  CoarseParameterRange = [10.0^k for k = -4:0.5:5]
  CoarseResultsBE = []
  for α in CoarseParameterRange
    Xstar, Iterates, FunIterates, TimeVector = OPT.tnSubgradientMethod(instanceBE, X0BE, Τ, α)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseSUBG-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  for α in CoarseParameterRange
    Xstar, Iterates, FunIterates, TimeVector = OPT.tnSubgradientMethod(instanceCA, X0CA, Τ, α)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseSUBG-CA.jld2", CoarseResultsCA)
end


function HPO_Coarse_SUBG_EP()
  CoarseParameterRange = [10.0^k for k = -4:0.5:5]
  CoarseResultsBE = []
  @info "Belgian"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tEstimatedPolyak(instanceBE, X0BE, Τ, α)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseSUBG-EP-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  @info "Californian"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tEstimatedPolyak(instanceCA, X0CA, Τ, α)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseSUBG-EP-CA.jld2", CoarseResultsCA)
end

function HPO_Coarse_LSG()
  CoarseParameterRange = [10.0^k for k = -4:0.5:5]
  CoarseResultsBE = []
  @info "Belgian"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tlastSubgradientMethod(instanceBE, X0BE, 900, α)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseLSG-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  @info "Californian"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tlastSubgradientMethod(instanceCA, X0CA, 600, α)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseLSG-CA.jld2", CoarseResultsCA)
end

function HPO_Coarse_BLM()
  CoarseParameterRange = [k for k in LinRange(0.2, 0.98, 10)]
  CoarseResultsBE = []
  @info "Running BE8"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBundleLevelMethod(instanceBE, X0BE, Τ, α)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//nCoarseBLM-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  @info "Running CA1"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBundleLevelMethod(instanceCA, X0CA, Τ, α)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//nCoarseBLM-CA.jld2", CoarseResultsCA)
end

function HPO_Coarse_BPLM_B()
  CoarseParameterRange = [k for k in LinRange(0.01, 0.95, 20)]
  CoarseResultsBE = []
  @info "Running BE8"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBundleProximalLevelMethod(instanceBE, X0BE, Τ, α)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseBPLM-B-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  @info "Running CA1"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBundleProximalLevelMethod(instanceCA, X0CA, Τ, α)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseBPLM-B-CA.jld2", CoarseResultsCA)
end

function HPO_Coarse_BPLM_L()
  CoarseParameterRange = [k for k in LinRange(0.01, 0.95, 20)]
  CoarseResultsBE = []
  @info "Running BE8"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBPLM(instanceBE, X0BE, Τ, α)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseBPLM-L-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  @info "Running CA1"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBPLM(instanceCA, X0CA, Τ, α)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseBPLM-L-CA.jld2", CoarseResultsCA)
end


function HPO_Coarse_FGM()
  CoarseParameterRange = [1e-6, 1e-5, 1e-4, 1e-3]
  CoarseResultsBE = []
  @info "Running BE8"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tShiftedFGM(instanceBE, X0BE, Τ, α, 1e-8)
    push!(CoarseResultsBE, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseFGM-BE.jld2", CoarseResultsBE)

  CoarseResultsCA = []
  @info "Running CA1"
  for α in CoarseParameterRange
    @info "Testing α = $α"
    Xstar, Iterates, FunIterates, TimeVector = OPT.tShiftedFGM(instanceCA, X0CA, Τ, α, 1e-8)
    push!(CoarseResultsCA, [α, Xstar, Iterates, FunIterates, TimeVector])
  end
  save_object("results//HPopt//CoarseFGM-CA.jld2", CoarseResultsCA)
end