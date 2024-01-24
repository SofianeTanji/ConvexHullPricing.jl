@info "Imports ..."
using Revise
using ConvexHullPricing
using DataFrames
using SparseArrays
using Plots
using JLD2
using LinearAlgebra, LaTeXStrings
using ProgressBars
const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer
BEinstances = []
for file in readdir(".//data//belgian"; join=true)
    push!(BEinstances, UT.load_data(file))
end
CAinstances = []
for file in readdir(".//data//ca"; join=true)
    push!(CAinstances, UT.load_data(file))
end

function make_it_monotone(array)
    new_array = Float64[]
    push!(new_array, array[1])
    for elt in array[2:end]
        if elt < minimum(new_array)
            push!(new_array, elt)
        else
            push!(new_array, minimum(new_array))
        end
    end
    return new_array
end


N_MIN = 15
τ = N_MIN * 60
idx = 1
@info "Instance $idx"
Belgian = BEinstances[idx]
X0 = UT.LP_Relaxation(Belgian)
FS = maximum(load_object("results//december//UltraOptimalRunsBE.jld2")[idx][3])
# Compilation
_, _, Hval, _, UB, LB = OPT.tBundleLevelMethod(Belgian, X0, 1., 0.5)
_, _, Hval, _, UB, LB = OPT.tBPLM(Belgian, X0, 1., 0.5)
_, _, Hval, _, UB, LB = OPT.tBundleProximalLevelMethod(Belgian, X0, 1., 0.5)

RanLambda = [k for k in LinRange(.01, .99, 99)]
BLM_vals = Float64[]
BPLM_vals = Float64[]
B_PLM_vals = Float64[]

BLM_UB = []
BPLM_UB = []
B_PLM_UB = []

BLM_LB = []
BPLM_LB = []
B_PLM_LB = []

cBPLM_vals = Float64[]
cBPLM_UB = []
cBPLM_LB = []


for λ in RanLambda
	@info "λ = $λ ; BLM run -->"
	_, _, Hval, _, UB, LB = OPT.tBundleLevelMethod(Belgian, X0, τ, λ)
	push!(BLM_vals, minimum(FS .- Hval))
	push!(BLM_UB, UB)
	push!(BLM_LB, LB)
end
save_object("BLM_tuning_on_Autumn_WD", [RanLambda, BLM_vals, BLM_UB, BLM_LB])
for λ in RanLambda
	@info "λ = $λ ; B(P)LM run -->"
	_, _, Hval, _, UB, LB = OPT.tBPLM(Belgian, X0, τ, λ)
	push!(B_PLM_vals, minimum(FS .- Hval))
	push!(B_PLM_UB, UB)
	push!(B_PLM_LB, LB)
end
save_object("B(P)LM_tuning_on_Autumn_WD", [RanLambda, B_PLM_vals, B_PLM_UB, B_PLM_LB])

for λ in RanLambda
	@info "λ = $λ ; BPLM run -->"
	_, _, Hval, _, UB, LB = OPT.tBundleProximalLevelMethod(Belgian, X0, τ, λ)
	push!(BPLM_vals, minimum(FS .- Hval))
	push!(BPLM_UB, UB)
	push!(BPLM_LB, LB)
end
save_object("BPLM_tuning_on_Autumn_WD", [RanLambda, BPLM_vals, BPLM_UB, BPLM_LB])


# I MAY HAVE CORRECTED BPLM

Ran2Lambda = [k for k in LinRange(.01, .99, 25)]
for λ in Ran2Lambda
	@info "λ = $λ ; BPLM run -->"
	_, _, Hval, _, UB, LB = OPT.tCorrectBundleProximalLevelMethod(Belgian, X0, τ, λ)
	push!(cBPLM_vals, minimum(FS .- Hval))
	push!(cBPLM_UB, UB)
	push!(cBPLM_LB, LB)
end
save_object("cBPLM_tuning_on_Autumn_WD", [Ran2Lambda, cBPLM_vals, cBPLM_UB, cBPLM_LB])
