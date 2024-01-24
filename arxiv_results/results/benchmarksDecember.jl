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

N_MIN = 0.1
τ = N_MIN * 60
idx = 1
@info "Instance $idx"
Bel = BEinstances[idx]
LPBE = UT.LP_Relaxation(Bel)
OptimalBE = maximum(load_object("results//december//UltraOptimalRunsBE.jld2")[idx][3])

_, Hit, Hval, _ = OPT.tShiftedFGM(Bel, LPBE, τ, 0.002, 1e-8)

N_MIN = 15
τ = N_MIN * 60
# GM = 0.003 FGM : 0.002
# R1 = 357 ; R2 = 21 ; R3 = 68.6 ; R4 = 13.4 ; R5 = 26.4 ; R6 = 18.6 ; R7 = 226.9 ; R8 = 121# 

for I = 1:8
    
end

for I=1:8
    @info "Instance $I"
    Belgian = BEinstances[I]
    X0 = UT.LP_Relaxation(Belgian)
    OptimalBE = maximum(load_object("results//december//UltraOptimalRunsBE.jld2")[I][3])
    @info "Ready to run FGM $I"
    _, Hit, Hval, Htime = OPT.tShiftedFGM(Belgian, X0, τ, 0.002, 1e-8)
    Htrue = Float64[]
    @info "done"
    for price in Hit
        push!(Htrue, UT.exact_oracle(Belgian, price)[1])
    end
    @info "FGM run $I fully complete."
    save_object("belgian$(I)FGM", [Htime, OptimalBE .- Htrue])
end
#=
for I=1:8
    @info "Instance $I"
    Belgian = BEinstances[I]
    X0 = UT.LP_Relaxation(Belgian)
    OptimalBE = maximum(load_object("results//december//UltraOptimalRunsBE.jld2")[I][3])
	@info "Opt = $OptimalBE"
    @info "Ready to run Pol $I"
    _, Hiter, Hval, Htime = OPT.tPolyakMethod(Belgian, X0, τ, OptimalBE)
    @info "T Pol run $I fully complete."
	@info "GAP is $(OptimalBE .- maximum(Hval))"
    save_object("belgian$(I)Pol", [Htime, OptimalBE .- Hval])
end

for I=16:20
    @info "Instance $I"
    Cali = CAinstances[I]
    X0 = UT.LP_Relaxation(Cali)
    OptimalCA = maximum(load_object("results//december//UltraOptimalRunsCA.jld2")[I][3])
    @info "Ready to run BPLM $I"
    _, Hiter, Hval, Htime = OPT.tBPLM(Cali, X0, τ, .98)
    @info "BPLM run $I fully complete."
    save_object("californian$(I)BPLM", [Htime, OptimalBE .- Hval])
end =#