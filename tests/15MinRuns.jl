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
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    F_star = maximum(load_object("results//optimal_values//NewOptRunBE$(idx).jld2")[3]) # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tPolyakMethod(instance, X0, BUDGET, F_star)
    @info "Done. Saving ..."
    save_object("results//15min_runs//PolyakMethodBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    F_star = maximum(load_object("results//optimal_values//NewOptRunCA$(idx).jld2")[3]) # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tPolyakMethod(instance, X0, BUDGET, F_star)
    @info "Done. Saving ..."
    save_object("results//15min_runs//PolyakMethodCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
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
  for (idx, instance) in enumerate(BEinstances)
    @info "Instance BE #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.2 # Parameter
    Xstar, Iterates, FunIterates, TimeVector = OPT.tBundleLevelMethod(instance, X0, BUDGET, α)
    @info "Done. Saving ..."
    save_object("results//15min_runs//BundleLevelMethodBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
  for (idx, instance) in enumerate(CAinstances)
    @info "Instance CA #$idx - Running for 15 minutes."
    X0 = UT.LP_Relaxation(instance) # First iterate
    α = 0.3 # Parameter
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