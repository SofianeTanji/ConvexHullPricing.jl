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
for file in readdir("data\\belgian"; join=true)
    push!(BEinstances, UT.load_data(file))
end
CAinstances = []
for file in readdir("data\\ca"; join=true)
    push!(CAinstances, UT.load_data(file))
end
RTSinstances = []
for file in readdir("data\\rts_gmlc"; join=true)
    push!(RTSinstances, UT.load_data(file))
end
GMLC = RTSinstances[2]
LPRTS = UT.LP_Relaxation(GMLC)
ObjMRTS = UT.Matching(GMLC).Obj

Belgian = BEinstances[1]
Californian = CAinstances[1]

LPBE = UT.LP_Relaxation(Belgian)
ObjMBE = UT.Matching(Belgian).Obj
LPCA = UT.LP_Relaxation(Californian)
ObjMCA = UT.Matching(Californian).Obj

OptimalRTS = maximum(load_object("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\tests\\UltraOptimalRunsRTS.jld2")[1][3])
OptUpliftsRTS = ObjMRTS - OptimalRTS


OptimalCA = maximum(load_object("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\tests\\UltraOptimalRunsCA.jld2")[1][3])
OptUpliftsCA = ObjMCA - OptimalCA

BUDGET = 0.5 * 60

# Bundle Level Method
OPT.tBundleLevelMethod(GMLC, LPRTS, 1, .5, 1) # Compilation purposes

range = [.35]  # (.2 7856; .4 7432)

# RTS
time_vecs = []
values = []
for alpha in range
    xstar, iterates, fiterates, timevector = OPT.tBundleLevelMethod(GMLC, LPRTS, BUDGET, alpha, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMRTS .- fiterates)
end
for i=1:5
    @info minimum(values[i]) - OptUpliftsRTS
end
# BELGIAN
time_vecs = []
values = []
for alpha in range
    xstar, iterates, fiterates, timevector = OPT.tBundleLevelMethod(Belgian, LPBE, BUDGET, alpha, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMBE .- fiterates)
end

for i=1:6
    @info minimum(values[i]) - OptUpliftsBE
end

# CA
time_vecsCA = []
valuesCA = []
for alpha in range
    xstar, iterates, fiterates, timevector = OPT.tBundleLevelMethod(Californian, LPCA, BUDGET, alpha, 1)
    push!(time_vecsCA, timevector)
    push!(valuesCA, ObjMCA .- fiterates)
end

for i=1:5
    @info minimum(valuesCA[i]) - OptUpliftsCA
end
time_vecsCA[1]
valuesCA[1]
plot(time_vecsCA[1][2:end], valuesCA[1] .- OptUpliftsCA, label = "alpha = $(range[1])")
plot!(time_vecsCA[2][2:end], valuesCA[2] .- OptUpliftsCA, label = "alpha = $(range[2])")
plot!(time_vecsCA[3][2:end], valuesCA[3] .- OptUpliftsCA, label = "alpha = $(range[3])")
plot!(time_vecsCA[4][2:end], valuesCA[4] .- OptUpliftsCA, label = "alpha = $(range[4])")
plot!(time_vecsCA[5][2:end], valuesCA[5] .- OptUpliftsCA, label = "alpha = $(range[5])")
plot!(yscale =:log, xlabel = "Time in seconds", ylabel = "Distance to optimum")

## CONCLUSION BE = 0.95, CA = 0.98


# Bundle Proximal Level Method
OPT.tBundleProximalLevelMethod(GMLC, LPRTS, 1, .5, 1) # Compilation purposes
range = [.99, .999]

# RTS
time_vecs = []
values = []
for alpha in range
    xstar, iterates, fiterates, timevector = OPT.tBundleProximalLevelMethod(GMLC, LPRTS, BUDGET, alpha, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMRTS .- fiterates)
end
for i=1:2
    @info minimum(values[i]) - OptUpliftsRTS
end

# Californian
time_vecs = []
values = []
for alpha in LinRange(.999, 1., 10)
    xstar, iterates, fiterates, timevector = OPT.tBundleProximalLevelMethod(Californian, LPCA, BUDGET, alpha, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMCA .- fiterates)
    @info "level is $alpha, distance to f* is: $(minimum((ObjMCA .- fiterates)) - OptUpliftsCA)"
end
plot(time_vecs[end][2:end], values[end] .- OptUpliftsCA, label = "Alpha = $(LinRange(.999, 1., 10)[end])")

plot(title="Hyperparameter Californian FGM")
for i=1:length(time_vecs)
    plot!(time_vecs[i][2:end], values[i] .- OptUpliftsCA, label = "Alpha = $(LinRange(.999, 1., 10)[i])")
end
plot!(xlabel = "Time seconds")

for i=1:20
    @info minimum(values[i]) - OptUpliftsCA
end

xstar, iterates, fiterates, timevector = OPT.tBundleProximalLevelMethod(Belgian, LPBE, 1200, .9999, 1)
_, _, vals, timevec = OPT.tBundleProximalLevelMethod(Belgian, LPBE, 300, .9999, 1)

plot(timevector[2:end], ObjMBE .- fiterates .- OptUpliftsBE, label = "pblm")
plot!(timevec[2:end], ObjMBE .- vals .- OptUpliftsBE, label = "Bigg alpha")
_, _, blmf, blmtime = OPT.tBundleLevelMethod(Belgian, LPBE, 1200, .999, 1)
plot!(blmtime[2:end], ObjMBE .- blmf .- OptUpliftsBE, label = "BLM", ylims = (0, 7e4))
plot(time_vecs[1][2:end], values[1] .- OptUpliftsBE, label = "alpha = $(range[1])")
plot!(time_vecs[2][2:end], values[2] .- OptUpliftsBE, label = "alpha = $(range[2])")
plot!(time_vecs[3][2:end], values[3] .- OptUpliftsBE, label = "alpha = $(range[3])")
plot!(time_vecs[4][2:end], values[4] .- OptUpliftsBE, label = "alpha = $(range[4])")
plot!(time_vecs[5][2:end], values[5] .- OptUpliftsBE, label = "alpha = $(range[5])")
plot!(time_vecs[6][2:end], values[6] .- OptUpliftsBE, label = "alpha = $(range[6])")

hline!([7566.31], label = "reached by blm in same time", ylims=(-1000,100000))
minimum(ObjMBE .- blmf .- OptUpliftsBE)
#plot!(time_vecs[5][2:end], values[5] .- OptUpliftsBE, label = "alpha = $(range[5])")
#plot!(time_vecs[6][2:end], values[6] .- OptUpliftsBE, label = "alpha = $(range[6])")

plot!(yscale =:log, xlabel = "Time in seconds", ylabel = "Distance to optimum")

# D-Adaptation et DowG
OPT.tDAdaptation(GMLC, LPRTS, 1, 1., 1e-7) # Compilation purposes
rangeD = [2e-8]

BUDGET = 5 * 60

# RTS
time_vecs = []
values = []
for alpha in rangeD
    xstar, iterates, fiterates, timevector = OPT.tDAdaptation(GMLC, LPRTS, BUDGET, alpha, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMRTS .- fiterates)
end
for i=1:1
    @info minimum(values[i]) - OptUpliftsRTS
end
idx = 1

plot(title="test")
for idx=1:1
    plot!(time_vecs[idx][2:end], values[idx], label = "$idx")
end
plot!(xlabel = "time")
# Belgian
time_vecs = []
values = []
for D in rangeD
    xstar, iterates, fiterates, timevector = OPT.tDAdaptation(Belgian, LPBE, 60, D, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMBE .- fiterates)
    @info "size is $D"
    @info minimum((ObjMBE .- fiterates)) - OptUpliftsBE
end

plot(time_vecs[1][2:end], values[1] .- OptUpliftsBE, label = "D = $(rangeD[1])")
plot(time_vecs[2][2:end], values[2] .- OptUpliftsBE, label = "D = $(rangeD[2])")
plot!(time_vecs[3][2:end], values[3] .- OptUpliftsBE, label = "D = $(rangeD[3])")
plot!(time_vecs[4][2:end], values[4] .- OptUpliftsBE, label = "D = $(rangeD[4])")
plot!(time_vecs[5][2:end], values[5] .- OptUpliftsBE, label = "D = $(rangeD[5])")
plot!(time_vecs[6][2:end], values[6] .- OptUpliftsBE, label = "D = $(rangeD[6])")
plot!(time_vecs[7][2:end], values[7] .- OptUpliftsBE, label = "D = $(rangeD[7])")


rangeD = [0.05, 0.1, 0.5]
# Californian
time_vecs = []
values = []
for D in rangeD
    xstar, iterates, fiterates, timevector = OPT.tDAdaptation(Californian, LPCA, 60, D, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMCA .- fiterates)
    @info "size is $D"
    @info minimum((ObjMCA .- fiterates)) - OptUpliftsCA
end

plot(time_vecs[1][2:end], values[1] .- OptUpliftsCA, label = "D = $(rangeD[1])")
plot!(time_vecs[2][2:end], values[2] .- OptUpliftsCA, label = "D = $(rangeD[2])")
plot!(time_vecs[3][2:end], values[3] .- OptUpliftsCA, label = "D = $(rangeD[3])")

# D-Adapt Belgian D = 20, Californian = 0.1

OPT.tDowG(GMLC, LPRTS, 1, 1., 1) # Compilation purposes

rangeD = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2]

#RTS
time_vecs = []
values = []
for D in rangeD
    xstar, iterates, fiterates, timevector = OPT.tDowG(GMLC, LPRTS, 60, D, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMRTS .- fiterates)
    @info "size is $D"
    @info minimum((ObjMRTS .- fiterates)) - OptUpliftsRTS
end
plot(title="test")
for idx=1:3
    plot!(time_vecs[idx][2:end], values[idx], label = "$idx")
end
plot!(xlabel = "time")
# Belgian
time_vecs = []
values = []
for D in rangeD
    xstar, iterates, fiterates, timevector = OPT.tDowG(Belgian, LPBE, 60, D, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMBE .- fiterates)
    @info "size is $D"
    @info minimum((ObjMBE .- fiterates)) - OptUpliftsBE
end
plot(time_vecs[1][2:end], values[1] .- OptUpliftsBE, label = "D = $(rangeD[1])")
plot!(time_vecs[2][2:end], values[2] .- OptUpliftsBE, label = "D = $(rangeD[2])")

rangeD = [1e-2, 1e-1]

# Californian
time_vecs = []
values = []
for D in rangeD
    xstar, iterates, fiterates, timevector = OPT.tDowG(Californian, LPCA, 60, D, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMCA .- fiterates)
    @info "size is $D"
    @info minimum((ObjMCA .- fiterates)) - OptUpliftsCA
end
plot(time_vecs[1][2:end], values[1] .- OptUpliftsCA, label = "D = $(rangeD[1])")
plot!(time_vecs[2][2:end], values[2] .- OptUpliftsCA, label = "D = $(rangeD[2])")
plot!(time_vecs[3][2:end], values[3] .- OptUpliftsCA, label = "D = $(rangeD[3])")
plot!(time_vecs[4][2:end], values[4] .- OptUpliftsCA, label = "D = $(rangeD[4])")

# DowG BE: D = 80, CA: D = 0.1

## SUBGRADIENT METHOD 

OPT.tSubgradientMethod(Belgian, LPBE, 1, 1., 1) # Compilation purposes

rangeG = [0.001, 0.005, 0.01, 0.05]
# Belgian
time_vecs = []
values = []
for G in rangeG
    xstar, iterates, fiterates, timevector = OPT.tSubgradientMethod(Belgian, LPBE, 60, G, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMBE .- fiterates)
    @info "size is $G"
    @info minimum((ObjMBE .- fiterates)) - OptUpliftsBE
end
plot(time_vecs[1][2:end], values[1] .- OptUpliftsBE, label = "D = $(rangeG[1])")
plot!(time_vecs[2][2:end], values[2] .- OptUpliftsBE, label = "D = $(rangeG[2])")
plot!(time_vecs[3][2:end], values[3] .- OptUpliftsBE, label = "D = $(rangeG[3])")
plot!(time_vecs[4][2:end], values[4] .- OptUpliftsBE, label = "D = $(rangeG[4])")
plot!(time_vecs[5][2:end], values[5] .- OptUpliftsBE, label = "D = $(rangeG[5])")
plot!(time_vecs[6][2:end], values[6] .- OptUpliftsBE, label = "D = $(rangeG[6])")

rangeD = [1e-2, 1e-1]

# Californian
time_vecs = []
values = []

rangeG = [0.00001, 0.0001, 0.001]
for G in rangeG
    xstar, iterates, fiterates, timevector = OPT.tSubgradientMethod(Californian, LPCA, 60, G, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMCA .- fiterates)
    @info "size is $G"
    @info minimum((ObjMCA .- fiterates)) - OptUpliftsCA
end
plot(time_vecs[1][2:end], values[1] .- OptUpliftsCA, label = "D = $(rangeG[1])")
plot(time_vecs[2][2:end], values[2] .- OptUpliftsCA, label = "D = $(rangeG[2])")
plot!(time_vecs[3][2:end], values[3] .- OptUpliftsCA, label = "D = $(rangeG[3])")
plot!(time_vecs[4][2:end], values[4] .- OptUpliftsCA, label = "D = $(rangeG[4])")

# SUBG : Belgian Gamma = 0.01 ; Californian Gamma = 0.0001

# ESTIMATED POLYAK

OPT.tEstimatedPolyak(Belgian, LPBE, 1, 1., 1) # Compilation purposes

rangeG = [5e5]
# Belgian
time_vecs = []
values = []
for G in rangeG
    xstar, iterates, fiterates, timevector = OPT.tEstimatedPolyak(Belgian, LPBE, 180, G, 1)
    push!(time_vecs, timevector)
    push!(values, ObjMBE .- fiterates)
    @info "size is $G"
    @info minimum((ObjMBE .- fiterates)) - OptUpliftsBE
end
plot(time_vecs[1][2:end], values[1] .- OptUpliftsBE, label = "D = $(rangeG[1])")
plot!(time_vecs[2][2:end], values[2] .- OptUpliftsBE, label = "D = $(rangeG[2])")
plot!(time_vecs[3][2:end], values[3] .- OptUpliftsBE, label = "D = $(rangeG[3])")
plot!(time_vecs[4][2:end], values[4] .- OptUpliftsBE, label = "D = $(rangeG[4])")
plot!(ylims=(0,78000))

rangeG = [200, 300, 400, 500, 600]

# Californian
timevecs, values = [], []
for G in rangeG
    _, _, fval, timevec = OPT.tEstimatedPolyak(Californian, LPCA, 60, G, 1)
    push!(timevecs, timevec)
    push!(values, ObjMCA .- fval)
    @info "gamma = $G and best distance : $(minimum((ObjMCA .- fval)) - OptUpliftsCA)"
end
plot(title="Hyperparameter Californian Est Polyak")
for i=1:length(timevecs)
    plot!(timevecs[i][2:end], values[i] .- OptUpliftsCA, label = "Gamma = $(rangeG[i])")
end
plot!(ylims=(0, 200))

# EST POLYAK BELGIAN = 5e5, Californian = 400

## STOCHASTIC GRADIENT METHOD
rangeEps = [1e-6, 1e-5, 1e-4]
_, its, _, timevec = OPT.tFastGradientMethod(Californian, LPCA, 10, 1e-6, 1)
fval = Float64[]
for price in its
    push!(fval, UT.exact_oracle(Californian, price)[1])
end
fval
_, _, myfs, _ = OPT.tEstimatedPolyak(Californian, LPCA, 10, 400, 1)
myfs
# Californian
timevecs, values = [], []
for eps in rangeEps
    _, its, _, timevec = OPT.tFastGradientMethod(Californian, LPCA, 120, eps, 1)
    push!(timevecs, timevec)
    fval = Float64[]
    for price in its
        push!(fval, UT.exact_oracle(Californian, price)[1])
    end
    push!(values, ObjMCA .- fval)
    @info "smoothing parameter = $eps and best distance : $(minimum((ObjMCA .- fval)) - OptUpliftsCA)"
end

plot(title="Hyperparameter Californian FGM")
for i=1:length(timevecs)
    plot!(timevecs[i][1:end], values[i] .- OptUpliftsCA, label = "Gamma = $(rangeEps[i])")
end
plot!(ylims=(0, 1e4))
