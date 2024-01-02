@info "Imports ..."
using Revise
using ConvexHullPricing
using DataFrames
using Plots
using JLD2
using ProgressBars
const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer
BEinstances = []
for file in readdir(".\\data\\belgian"; join=true)
    push!(BEinstances, UT.load_data(file))
end
CAinstances = []
for file in readdir(".\\data\\ca"; join=true)
    push!(CAinstances, UT.load_data(file))
end

N_MIN = 1
BUDGET = N_MIN * 30

Belgian = BEinstances[1]
LPBE = UT.LP_Relaxation(Belgian)
_, SUBGiter, SUBGval, SUBGtime = OPT.tSubgradientMethod(Belgian, LPBE, BUDGET, 0.01)
_, H3iter, H3val, H3time = OPT.tRSG(Belgian, LPBE, BUDGET, 3, 0.025, -1)
_, Piter, Pval, Ptime = OPT.tPolyakMethod(Belgian, LPBE, BUDGET, OptimalBE, 1)

make_it_monotone(OptimalBE .- Pval)
make_it_monotone(OptimalBE .- SUBGval)
SUBGval

plot(Ptime, OptimalBE .- Pval, label = "neg Polyak")
plot!(nPtime, OptimalBE .- nPval, label = "pos Polyak")

plot!(H3time[2:37], OptimalBE .- H3val[1:36], label = "RSG")
H3time[2:37]
H3val[1:37]
for I=1:8
    @info "Instance $I"
    Belgian = BEinstances[I]
    LPBE = UT.LP_Relaxation(Belgian)
    OptimalBE = maximum(load_object("results\\december\\UltraOptimalRunsBE.jld2")[I][3])
    @info "Ready to run RSG $I"
    _, Hiter, Hval, Htime = OPT.tRSG(Belgian, LPBE, BUDGET, 3, 0.025, 1)
    @info "RSG run $I fully complete."
    save_object("belgian$(I)RSG", [Htime, OptimalBE .- Hval])
end

for I=1:8
    @info "Instance $I"
    Belgian = BEinstances[I]
    LPBE = UT.LP_Relaxation(Belgian)
    OptimalBE = maximum(load_object("results\\december\\UltraOptimalRunsBE.jld2")[I][3])
    @info "Ready to run Polyak $I"
    _, Hiter, Hval, Htime = OPT.tPolyakMethod(Belgian, LPBE, BUDGET, OptimalBE)
    @info "POL run $I fully complete."
    save_object("belgian$(I)Pol", [Htime, OptimalBE .- Hval])
end


for I=1:8
    @info "Instance $I"
    Belgian = BEinstances[I]
    LPBE = UT.LP_Relaxation(Belgian)
    OptimalBE = maximum(load_object("results\\december\\UltraOptimalRunsBE.jld2")[I][3])
    @info "Ready to run CG $I"
    _, CGiter, CGtime = OPT.tColumnGeneration(Belgian, LPBE, BUDGET, 1e-5)
    @info "CG run $I complete. Evaluating dual values"
    CGval = Float64[]
    for elt in CGiter
        push!(CGval, UT.exact_oracle(Belgian, elt)[1])
    end
    @info "CG run $I fully complete."
    save_object("belgian$(I)CG", [CGtime, OptimalBE .- CGval])
end
for I=16:20
    @info "Instance $I"
    Belgian = CAinstances[I]
    LPBE = UT.LP_Relaxation(Belgian)
    OptimalBE = maximum(load_object("results\\december\\UltraOptimalRunsCA.jld2")[I][3])
    @info "Ready to run CG $I"
    _, CGiter, CGtime = OPT.tColumnGeneration(Belgian, LPBE, BUDGET, 1e-5)
    @info "CG run $I complete. Evaluating dual values"
    CGval = Float64[]
    for elt in CGiter
        push!(CGval, UT.exact_oracle(Belgian, elt)[1])
    end
    @info "CG run $I fully complete."
    save_object("californian$(I)CG", [CGtime, OptimalBE .- CGval])
end

"""
N_MIN = 2
BUDGET = N_MIN * 60

_, BLMiter, BLMval, BLMtime = OPT.tBundleLevelMethod(Belgian, LPBE, BUDGET, 0.35, 1)
_, BPLMiter, BPLMval, BPLMtime = OPT.tBundleProximalLevelMethod(Belgian, LPBE, BUDGET, .99, 1)
_, DAiter, DAval, DAtime = OPT.tDAdaptation(Belgian, LPBE, BUDGET, .1)
_, DOWGiter, DOWGval, DOWGtime = OPT.tDowG(Belgian, LPBE, BUDGET, .1)
_, ESTPOLiter, ESTPOLval, ESTPOLtime = OPT.tEstimatedPolyak(Belgian, LPBE, BUDGET, 400)
_, SUBGiter, SUBGval, SUBGtime = OPT.tSubgradientMethod(Belgian, LPBE, BUDGET, 0.0001)

plot(title = "Comparison of methods for computing CH prices", xlabel = "Time in seconds", ylabel = "Function gap value", ylims = (0, 1e4))
plot!(BLMtime[2:end], OptimalBE .- BLMval, label = "BLM")
plot!(BPLMtime[2:end], OptimalBE .- BPLMval, label = "BPLM")
plot!(DAtime[2:end], OptimalBE .- DAval, label = "DA")
plot!(DOWGtime[2:end], OptimalBE .- DOWGval, label = "DOWG")
plot!(ESTPOLtime[2:end], OptimalBE .- ESTPOLval, label = "Est. Polyak")
plot!(SUBGtime[2:end], OptimalBE .- SUBGval, label = "SUBG")
savefig("belgian$I")
save_object("belgian$(I)CG", [BLMtime, OptimalBE .- BLMval])

save_object("belgian$(I)BLM", [BLMtime, OptimalBE .- BLMval])
save_object("belgian$(I)BPLM", [BPLMtime, OptimalBE .- BPLMval])
save_object("belgian$(I)DA", [DAtime, OptimalBE .- DAval])
save_object("belgian$(I)DOWG", [DOWGtime, OptimalBE .- DOWGval])
save_object("belgian$(I)ESTPOL", [ESTPOLtime, OptimalBE .- ESTPOLval])
save_object("belgian$(I)SUBG", [SUBGtime, OptimalBE .- SUBGval])
"""