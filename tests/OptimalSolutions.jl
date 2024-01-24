using Revise
using ConvexHullPricing
using DataFrames
using Plots
using JLD2
using ProgressBars, JuMP, Gurobi

const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer

#= BEinstances = []
for file in readdir("data//belgian"; join=true)
  push!(BEinstances, UT.load_data(file))
end

@info "Large and long computing ready to run !"

for (idx, instance) in enumerate(BEinstances)
    @info "Instance #$idx"
    X0 = UT.LP_Relaxation(instance)
    Xstar, Iterates, FunIterates, TimeVector = OPT.tOptimal(instance, X0, .4)
    save_object("NewOptRunBE$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
end 

CAinstances = []
for file in readdir("data//ca"; join=true)
  push!(CAinstances, UT.load_data(file))
end

for (idx, instance) in enumerate(CAinstances)
  @info "Instance #$idx"
  X0 = UT.LP_Relaxation(instance)
  Xstar, Iterates, FunIterates, TimeVector = OPT.tOptimal(instance, X0, 0.8)
  save_object("NewOptRunCA$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
end
=#
function OptimalRTS()
  RTS_GMLCinstances = []
  for file in readdir("data//rts_gmlc"; join=true)
    push!(RTS_GMLCinstances, UT.load_rts_data(file))
  end

  for (idx, instance) in enumerate(RTS_GMLCinstances)
    @info "Instance #$idx"
    X0 = UT.LP_Relaxation(instance)
    Xstar, Iterates, FunIterates, TimeVector = OPT.tOptimal(instance, X0, 0.4)
    save_object("NewOptRunRTSGMLC$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end

function OptimalFERC()
  FERCinstances = []
  for file in readdir("data//ferc"; join=true)
    push!(FERCinstances, UT.load_ferc_data(file))
  end

  for (idx, instance) in enumerate(FERCinstances)
    @info "Instance #$idx"
    X0 = UT.LP_Relaxation(instance)
    Xstar, Iterates, FunIterates, TimeVector = OPT.tOptimal(instance, X0, 0.4)
    save_object("NewOptRunFERC$(idx).jld2", [Xstar, Iterates, FunIterates, TimeVector])
  end
end