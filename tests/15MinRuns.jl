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
BUDGET = 15 * 60
@info "Loading done."
# METHOD 1 : Polyak Method
function RunPolyak()
  @info "Polyak Method"
  _, _, _, _ = OPT.tPolyakMethod(BEinstances[1], zeros(96), 0.01, 1.)
  @info "Precompilation done."
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    F_star = maximum(load_object("results//optimal_values//NewRefinedOptRunBE$(idx).jld2")[3]) # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tPolyakMethod(instance, X0, BUDGET, F_star)
    @info "Done. Saving ..."
    save_object("results//15min_runs//nPolyakMethodBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    F_star = maximum(load_object("results//optimal_values//NewRefinedOptRunCA$(idx).jld2")[3]) # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tPolyakMethod(instance, X0, BUDGET, F_star)
    @info "Done. Saving ..."
    save_object("results//15min_runs//nPolyakMethodCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end
# METHOD 2 : D-Adaptation
function RunDAdaptation()
  @info "D-Adaptation"
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 10. # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDAdaptation(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//D-AdaptationBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.15 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDAdaptation(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//D-AdaptationCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

# METHOD 3 : BLM
function RunBLM()
  @info "Bundle Level Method"
  Xstar, Iterates, FunIterates, TimeVector = OPT.tBundleLevelMethod(BEinstances[1], zeros(96), 1., 0.9)
  @info "Compiled. Now running."
  for (idx, instance) in enumerate(BEinstances)
      @info "Instance BE #$idx - Running for 15 minutes."
      X0 = UT.LP_Relaxation(instance) # First iterate
      α = 0.9 # Parameter
      Xstar, Iterates, FunIterates, TimeVector = OPT.tBundleLevelMethod(instance, X0, BUDGET, α)
      @info "Done. Saving ..."
      save_object("results//15min_runs//BundleLevelMethodBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.9 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBundleLevelMethod(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//BundleLevelMethodCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

# METHOD 3 : BLM
function RunDowG()
  @info "DowG"
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 20. # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDowG(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//DowGBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.1 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tDowG(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//DowGCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

function RunSUBG()
  @info "SUBG"
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 100.0 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tnSubgradientMethod(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//SubGBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.3 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tnSubgradientMethod(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//SubGCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

function RunFGM()
  @info "FGM"
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 1e-4 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tShiftedFGM(instance, X0, BUDGET, α, 1e-8)
    @info "Done. Saving ..."
    save_object("results//15min_runs//FGMBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
      @info "Instance CA #$idx - Running for 15 minutes."
      X0 = UT.LP_Relaxation(instance) # First iterate
      α = 1e-5 # Parameter
      Xstar, Iterates, FunIterates, TimeVector = OPT.tShiftedFGM(instance, X0, BUDGET, α, 1e-8)
      @info "Done. Saving ..."
      save_object("results//15min_runs//FGMCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

function RunESTPOL()
  @info "EST POL"
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 32000. # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tEstimatedPolyak(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//SubG-EP-BE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 300. # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tEstimatedPolyak(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//SubG-EP-CA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

# METHOD 3 : BLM
function RunBPLM_L()
  @info "Bundle Proximal Level Method - Last Iterate"
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.95 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBPLM(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//BundleProximalLevelMethod-L-BE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.9 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBPLM(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//BundleProximalLevelMethod-L-CA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

function RunDW()
  @info "Column Generation"
  Xstar, Iterates, FunIterates, TimeVector = OPT.tColumnGeneration(BEinstances[1], zeros(96), 5., 1e-7)
  @info "Compiled. Now starting."
  verb = 1
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 1e-5 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tColumnGeneration(instance, X0, BUDGET, α, verb)
    @info "Done. Saving ..."
    save_object("results//15min_runs//ColumnGeneration-BE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 1e-5 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tColumnGeneration(instance, X0, BUDGET, α, verb)
    @info "Done. Saving ..."
    save_object("results//15min_runs//ColumnGeneration-CA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

function RunCRG()
  @info "Column Row Generation"
  Xstar, Iterates, FunIterates, TimeVector = OPT.tCRG(BEinstances[1], 5., 1e-5)
  @info "Compiled. Now starting."
  verb = 1
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    # X0 = UT.LP_Relaxation(instance) # First iterate
    α = 1e-5 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tCRG(instance, BUDGET, α, verb)
    @info "Done. Saving ..."
    save_object("results//15min_runs//CRG-BE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    # X0 = UT.LP_Relaxation(instance) # First iterate
    α = 1e-5 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tCRG(instance, BUDGET, α, verb)
    @info "Done. Saving ..."
    save_object("results//15min_runs//CRG-CA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

function RunLSUBG()
  @info "Last-iterate optimal subgradient method"
  Xstar, Iterates, FunIterates, TimeVector = OPT.tlastSubgradientMethod(BEinstances[1], zeros(96), 1, 0.9)
  @info "Compiled. Now running."
  for (idx, instance) in enumerate(BEinstances)
      @info "Instance BE #$idx - Running for 15 minutes."
      X0 = UT.LP_Relaxation(instance) # First iterate
      α = 40. # Parameter
      Xstar, Iterates, FunIterates, TimeVector = OPT.tlastSubgradientMethod(instance, X0, 1200, α)
      @info "Done. Saving ..."
      save_object("results//15min_runs//LSUBGBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.2 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tlastSubgradientMethod(instance, X0, 400, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//LSUBGCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end