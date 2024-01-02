using Revise
using ConvexHullPricing
using DataFrames
using Plots
using JLD2
using ProgressBars, JuMP, Gurobi

const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer

RTSinstances = []
for file in readdir("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\data\\rts_gmlc"; join=true)
    push!(RTSinstances, UT.load_data(file))
end
UltraOptimalRuns = []
@info "Large and long computing ready to run !"

for (idx, instance) in enumerate(RTSinstances)
    @info "Instance #$idx"
    LP_Relax = UT.LP_Relaxation(instance)
    BLMxstar, BLMiterates, BLMfvalues, BLMtimevector = OPT.tOptimal(instance, LP_Relax, .2)
    push!(UltraOptimalRuns, [BLMxstar, BLMiterates, BLMfvalues, BLMtimevector])
end

RTSinstances
for i in [5, 8, 9, 10, 11]
    @info "Instance #$i"
    LPRTS = 40 * ones(48)
    BLMxstar, BLMiterates, BLMfvalues, BLMtimevector = OPT.tOptimal(RTSinstances[i], LPRTS, .35)
    UltraOptimalRuns[i] = [BLMxstar, BLMiterates, BLMfvalues, BLMtimevector]
end

idx = 12
plot(UltraOptimalRuns[idx][4][2:end], - UltraOptimalRuns[idx][3] .- minimum(.-UltraOptimalRuns[idx][3]), xlabel = "Time (s)", ylabel = "Distance to optimum", ylims=(0, 5e4), label = "RTS instance $idx", title="RTS dataset")

# 5, 8, 9, 10, 11, 14
save_object("UltraOptimalRunsRTS.jld2", UltraOptimalRuns)