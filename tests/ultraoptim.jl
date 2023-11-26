using Revise
using ConvexHullPricing
using DataFrames
using Plots
using JLD2
using ProgressBars, JuMP, Gurobi

const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer

BEinstances = []
for file in readdir("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\data\\belgian"; join=true)
    push!(BEinstances, UT.load_data(file))
end
CAinstances = []
for file in readdir("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\data\\ca"; join=true)
    push!(CAinstances, UT.load_data(file))
end
allInstances = CAinstances
UltraOptimalRuns = []
@info "Large and long computing ready to run !"
for instance in allInstances
    @info "New instance"
    LP_Relax = UT.LP_Relaxation(instance)
    BLMxstar, BLMiterates, BLMfvalues, BLMtimevector = OPT.tOptimal(instance, LP_Relax, .7)
    push!(UltraOptimalRuns, [BLMxstar, BLMiterates, BLMfvalues, BLMtimevector])
end

save_object("UltraOptimalRunsCA.jld2", UltraOptimalRuns)