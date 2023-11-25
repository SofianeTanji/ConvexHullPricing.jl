@info "Imports ..."
using Revise
using ConvexHullPricing
using DataFrames
using Plots
using JLD2
using ProgressBars

N_ITER = 1
BUDGET = 5 * 60
const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer
@info "Loading toy instance for precompilation"
ThermalGen = ConvexHullPricing.Utilitaries.ThermalGen(
    MinRunCapacity = [6],
    MaxRunCapacity = [16],
    RampUp = [5],
    RampDown = [5],
    StartUp = [6],
    ShutDown = [6],
    UpTime = [1],
    DownTime = [1],
    NoLoadConsumption = [10],
    MarginalCost = [53],
    FixedCost = [30],
)
instance = ConvexHullPricing.Utilitaries.Instance(
    LostLoad = 3000,
    Load = [6 11 16 11 11 16 11 16 20],
    ThermalGen = ThermalGen
)
LP_Relax = UT.LP_Relaxation(instance)
ObjM = UT.Matching(instance).Obj
ValueAtInitialPoint = UT.fast_oracle(instance, LP_Relax)[1] .- ObjM

@info UT.smooth_oracle(instance, LP_Relax, 1e-3)[1]
@info UT.super_fast_oracle(instance, LP_Relax)[1]

OPT.tBundleLevelMethod(instance, LP_Relax, N_ITER, .9)
OPT.tBundleProximalMethod(instance, LP_Relax, N_ITER, .9, 1e-1)
OPT.tDAdaptation(instance, LP_Relax, N_ITER)
OPT.tDowG(instance, LP_Relax, N_ITER)
OPT.tFastGradientMethod(instance, LP_Relax, N_ITER, 1e-5)
OPT.tEstimatedPolyak(instance, LP_Relax, N_ITER, 1.)
OPT.tSubgradientMethod(instance, LP_Relax, N_ITER, 1.)
OPT.tStochasticSubgradientMethod(instance, LP_Relax, N_ITER, 1., 1)
OPT.ColumnGeneration(instance, LP_Relax, N_ITER, 1e-5)
@info "All precompilation completed. Loading belgian instances"
BEinstances = []
for file in readdir("data\\belgian"; join=true)
    push!(BEinstances, UT.load_data(file))
end

function all_results4subplots(instances, budget)
    BLM_results = []
    BPM_results = []
    DA_results = []
    DowG_results = []
    ePol_results = []
    GD_results = []
    CG_results = []
    # CRG_results = []
    idx = 0
    for instance in instances
        idx += 1
        @info "[Instance #$idx]"
        LP_Relax = UT.LP_Relaxation(instance)
        ObjM = UT.Matching(instance).Obj
        @info "[BLM]"
        xs, it, fit, tvec = OPT.tBundleLevelMethod(instance, LP_Relax, budget, .99)
        push!(BLM_results, [ObjM .- fit, tvec[2:end]])
        @info "[BPM]"
        xs, it, fit, tvec = OPT.tBundleProximalMethod(instance, LP_Relax, budget, .98, 1e4)
        push!(BPM_results, [ObjM .- fit, tvec[2:end]])
        @info "[DA]"
        xs, it, fit, tvec = OPT.tDAdaptation(instance, LP_Relax, budget, 5e-2)
        push!(DA_results, [ObjM .- fit, tvec[2:end]])
        @info "[DowG]"
        xs, it, fit, tvec = OPT.tDowG(instance, LP_Relax, budget, 1e-4)
        push!(DowG_results, [ObjM .- fit, tvec[2:end]])
        @info "[FGM]"
        # xs, it, fit, tvec = OPT.tFastGradientMethod(instance, LP_Relax, budget, 1.)
        # push!(FGM_results, [ObjM .- fit, tvec[2:end]])
        @info "[ePol]"
        xs, it, fit, tvec = OPT.tEstimatedPolyak(instance, LP_Relax, budget, 1.)
        push!(ePol_results, [ObjM .- fit, tvec[2:end]])
        @info "[GD]"
        xs, it, fit, tvec = OPT.tSubgradientMethod(instance, LP_Relax, budget, 5e-4)
        push!(GD_results, [ObjM .- fit, tvec[2:end]])
        @info "[SGD]"
        # xs, it, fit, tvec = OPT.tStochasticSubgradientMethod(instance, LP_Relax, budget, 1e-3, 5)
        # push!(SGD_results, [ObjM .- fit, tvec[2:end]])
        @info "[CG]"
        xs, it, tvec = OPT.tColumnGeneration(instance, LP_Relax, budget, 1e-5)
        fit = Float64[]
        for price in it
            push!(fit, UT.super_fast_oracle(instance, price)[1])
        end
        push!(CG_results, [ObjM .- fit, tvec[2:end]])
        # @info "[CRG]"
        # xs, it, fit, tvec = OPT.tColumnAndRowGeneration(instance, budget, 1e-5)
        # push!(CRG_results, [ObjM .- fit, tvec[2:end]])
    end
    return BLM_results, BPM_results, DA_results, DowG_results, ePol_results, GD_results, CG_results#, CRG_results
end

BLM_results, BPM_results, DA_results, DowG_results, FGM_results, ePol_results, GD_results, SGD_results, CG_results = all_results4subplots(BEinstances, 300)

save_object("results//figures//beBLM.jld2", BLM_results)
save_object("results//figures//beBPM.jld2", BPM_results)
save_object("results//figures//beDA.jld2", DA_results)
save_object("results//figures//beDowG.jld2", DowG_results)
save_object("results//figures//beFGM.jld2", FGM_results)
save_object("results//figures//beePol.jld2", ePol_results)
save_object("results//figures//beGD.jld2", GD_results)
save_object("results//figures//beSGD.jld2", SGD_results)
save_object("results//figures//beCG.jld2", CG_results)

CAinstances = []
for file in readdir("data\\ca"; join=true)
    push!(CAinstances, UT.load_data(file))
end

CAinstances = CAinstances[11:20]

CAinstances

BLM_results, BPM_results, DA_results, DowG_results, ePol_results, GD_results, CG_results = all_results4subplots(CAinstances, 300)

save_object("results//figures//cabisBLM.jld2", BLM_results)
save_object("results//figures//cabisBPM.jld2", BPM_results)
save_object("results//figures//cabisDA.jld2", DA_results)
save_object("results//figures//cabisDowG.jld2", DowG_results)
save_object("results//figures//cabisePol.jld2", ePol_results)
save_object("results//figures//cabisGD.jld2", GD_results)
save_object("results//figures//cabisCG.jld2", CG_results)

plot(BLM_results[3][2], BLM_results[3][1], label = "BLM")
plot(BPM_results[3][2], BPM_results[3][1], label = "BPM")
plot!(DA_results[3][2], DA_results[3][1], label = "DA")
plot!(DowG_results[3][2], DowG_results[3][1], label = "DowG")
plot!(ePol_results[3][2], ePol_results[3][1], label = "ePol")
plot!(GD_results[3][2], GD_results[3][1], label = "GD")
plot!(CG_results[3][2], CG_results[3][1], label = "CG")