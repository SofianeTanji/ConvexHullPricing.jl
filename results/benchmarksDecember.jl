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
for file in readdir(".\\data\\belgian"; join=true)
    push!(BEinstances, UT.load_data(file))
end
CAinstances = []
for file in readdir(".\\data\\ca"; join=true)
    push!(CAinstances, UT.load_data(file))
end

function set_T(UT,T_max, NbGen)
	index_gen = zeros(NbGen)
	for g=1:NbGen
		index = 0
		# Intervals for generator "on" on prior, intervals [0,b]
		# index+=T_max+1; # [0,b] for b=0:T_max
		# Intervals [a,b] such that 1 <= a <= a+(UT-1) <= b <= T_max
		for a=1:T_max
			for b=a+UT[g]-1:T_max
				if(1 <= a && a <= a+UT[g]-1 && a+UT[g]-1 <= b && b <= T_max)
					index += 1
				end
			end
		end
		# Intervals for generator "on" past time T_max, intervals [a,T_max+1]
		# index+=T_max+1; # [a,T_max+1] for a=1:T_max+1
		# Interval for all period [0,T_max+1]
		# index+=1;
		index_gen[g] = index
	end
	A = sparse(zeros(NbGen, convert(Int64,maximum(index_gen))))
	B = sparse(zeros(NbGen, convert(Int64,maximum(index_gen))))

	for g=1:NbGen
		index = 1
		for a=1:T_max
			for b=a+UT[g]-1:T_max
				if(1 <= a && a <= a+UT[g]-1 && a+UT[g]-1 <= b && b <= T_max)
					A[g,index] = a
					B[g,index] = b
					index+=1
				end
			end
		end
	end
	A = convert(Matrix{Int64}, A)
	B = convert(Matrix{Int64}, B)
	index_gen = convert(Array{Int64}, index_gen)
	return (A,B,index_gen)
end

N_MIN = 15
τ = N_MIN * 60

idx = 1
@info "Instance $idx"
Bel = BEinstances[idx]
LPBE = UT.LP_Relaxation(Bel)
OptimalBE = maximum(load_object("results\\december\\UltraOptimalRunsBE.jld2")[idx][3])

BE1BLMtime, BE1BLMval = load_object("results\\december\\belgian1BLM")
BE1BPLMtime, BE1BPLMval = load_object("results\\december\\belgian1BPLM")
BE1CGtime, BE1CGval = load_object("results\\december\\belgian1CG")
BE1DAtime, BE1DAval = load_object("results\\december\\belgian1DA")
BE1DOWGtime, BE1DOWGval = load_object("results\\december\\belgian1DOWG")
BE1ESTPOLtime, BE1ESTPOLval = load_object("results\\december\\belgian1ESTPOL")

ArrOpt = load_object("results\\december\\UltraOptimalRunsBE.jld2")[idx]
R = norm(LPBE - ArrOpt[1])

# R1 = 357 ; R2 = 21 ; R3 = 68.6 ; R4 = 13.4 ; R5 = 26.4 ; R6 = 18.6 ; R7 = 226.9 ; R8 = 121# 
xngd, itngd, fngd, tngd = OPT.tnSubgradientMethod(Bel, LPBE, τ, 25)
length(itngd)
xnrgd, itnrgd, fnrgd, tnrgd = OPT.tnRSG(Bel, LPBE, τ, 3, 25)
xlngd, itlngd, flngd, tlngd = OPT.tlastSubgradientMethod(Bel, LPBE, length(itngd), R)
Smallestxlngd, Smallestitlngd, Smallestflngd, Smallesttlngd = OPT.tlastSubgradientMethod(Bel, LPBE, length(itngd), 13.4)
Meanxlngd, Meanitlngd, Meanflngd, Meantlngd = OPT.tlastSubgradientMethod(Bel, LPBE, length(itngd), 106.6)
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
plot(tlngd[2:end], make_it_monotone(OptimalBE .- flngd), label = "R for dataset", markershape = :circle, markersize = 0.1, ylims = (0, 8e4))
plot!(Meantlngd[2:end], make_it_monotone(OptimalBE .- Meanflngd), label = "mean R", markershape = :circle, markersize = 0.1, ylims = (0, 8e4))
plot!(BE1BLMtime[2:end], make_it_monotone(BE1BLMval), label = "Bundle Level Method", markershape = :circle, markersize = 0.1)
plot!(BE1BPLMtime[2:end], make_it_monotone(BE1BPLMval), label = "Bundle Proximal Level Method", markershape = :circle, markersize = 0.1)
plot!(BE1DAtime[2:end], make_it_monotone(BE1DAval), label = "D-Adaptation", markershape = :circle, markersize = 0.1)
plot!(BE1DOWGtime[2:end], make_it_monotone(BE1DOWGval), label = "DowG", markershape = :circle, markersize = 0.1)
plot!(BE1ESTPOLtime[2:end], make_it_monotone(BE1ESTPOLval), label = "Estimated Polyak", markershape = :circle, markersize = 0.1)
plot!(BE1CGtime[2:end], make_it_monotone(BE1CGval), label = "Dantzig-Wolfe Method", markershape = :circle, markersize = 0.1)
plot!(xlabel = "Time in seconds", ylabel = L"f(x_t) - f_*", title = "Comparison of different strategies for subgradient method - BE Autumn WD", margin = 1.1Plots.cm, size = (1400, 800))



@info "test"
plot(tlngd[2:end], OptimalBE .- flngd, label = "Last iterate Subgradient method", markershape = :circle, markersize = 0.1, xlims = (0,900), ylims = (0, 8e4))
plot!(tnrgd[2:end], OptimalBE .- fnrgd, label = "Restarted Subgradient method with normalized subgradients", markershape = :circle, markersize = 0.1)
plot!(tngd[2:end], OptimalBE .- fngd, label = "Subgradient method with normalized subgradients", markershape = :circle, markersize = 0.1)
plot!(BE1BLMtime[2:end], BE1BLMval, label = "Bundle Level Method", markershape = :circle, markersize = 0.1)

plot!(xlabel = "Time in seconds", ylabel = L"f(x_t) - f_*", title = "Comparison of different strategies for subgradient method - BE Autumn WD", margin = 1.1Plots.cm, size = (1400, 800))
savefig("MethodsBE1_try.pdf")

xgd, itgd, fgd, tgd = OPT.tSubgradientMethod(Bel, LPBE, τ, 0.01)
xrgd, itrgd, frgd, trgd = OPT.tRSG(Bel, LPBE, τ, 3, 0.01)

using LaTeXStrings
plot(tgd[2:end], OptimalBE .- fgd, label = "Subgradient method", markershape = :circle, markersize = 1)
plot!(tngd[2:end], OptimalBE .- fngd, label = "Subgradient method with normalized subgradients", markershape = :circle, markersize = 1)
plot!(trgd[2:end], OptimalBE .- frgd, label = "Restarted Subgradient method", markershape = :circle, markersize = 1)
plot!(tnrgd[2:end], OptimalBE .- fnrgd, label = "Restarted Subgradient method with normalized subgradients", markershape = :circle, markersize = 1)
plot!(xlabel = "Time in seconds", ylabel = L"f(x_t) - f_*", title = "Comparison of different strategies for subgradient method - BE Autumn WD", margin = 1.1Plots.cm, size = (1400, 800))
savefig("SubgradientMethods-autumnWD.pdf")


_, Hiter, Hval, Htime = OPT.tCRG(Bel, 20, 1e-5)

for I=1:8
    @info "Instance $I"
    Belgian = BEinstances[I]
    X0 = UT.LP_Relaxation(Belgian)
    OptimalBE = maximum(load_object("results\\december\\UltraOptimalRunsBE.jld2")[I][3])
    @info "Ready to run BPLM $I"
    _, Hiter, Hval, Htime = OPT.tBPLM(Belgian, X0, τ, .95)
    @info "BPLM run $I fully complete."
    save_object("belgian$(I)BPLM", [Htime, OptimalBE .- Hval])
end

for I=1:20
    @info "Instance $I"
    Cali = CAinstances[I]
    X0 = UT.LP_Relaxation(Cali)
    OptimalCA = maximum(load_object("results\\december\\UltraOptimalRunsCA.jld2")[I][3])
    @info "Ready to run BPLM $I"
    _, Hiter, Hval, Htime = OPT.tBPLM(Cali, X0, τ, .95)
    @info "BPLM run $I fully complete."
    save_object("californian$(I)BPLM", [Htime, OptimalBE .- Hval])
end