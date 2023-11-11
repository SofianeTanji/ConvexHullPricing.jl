using DataFrames, JLD2, ConvexHullPricing

beBLMdf = load_object("results//tables//beBLMdf.jld2")
caBLMdf = load_object("results//tables//caBLMdf.jld2")
BLMdf = vcat(beBLMdf, caBLMdf)

beBPMdf = load_object("results//tables//beBPMdf.jld2")
caBPMdf = load_object("results//tables//caBPMdf.jld2")
BPMdf = vcat(beBPMdf, caBPMdf)

beDAdf = load_object("results//tables//beDAdf.jld2")
caDAdf = load_object("results//tables//caDAdf.jld2")
DAdf = vcat(beDAdf, caDAdf)

beDowGdf = load_object("results//tables//beDowGdf.jld2")
caDowGdf = load_object("results//tables//caDowGdf.jld2")
DowGdf = vcat(beDowGdf, caDowGdf)

bePolyakdf = load_object("results//tables//bePolyakdf.jld2")
caPolyakdf = load_object("results//tables//caPolyakdf.jld2")
Polyakdf = vcat(bePolyakdf, caPolyakdf)

list_of_matchings = []
for instance in BLMdf.instance
    push!(list_of_matchings, ConvexHullPricing.Utilitaries.Matching(instance).Obj)
end

Polyakdf[!, :uplifts] .= Polyakdf[!, :fiterates]

for i=1:28
    Polyakdf.uplifts[i] = -(Polyakdf.uplifts[i] .- list_of_matchings[i])
end

for i=1:28
    @info "range of instance $i is [$(minimum(Polyakdf.uplifts[i])); $(maximum(Polyakdf.uplifts[i]))]"
end

save_object("results//objects//FinalBPMDF.jld2", BPMdf)

uplits = []
for arr in Polyakdf.uplifts
    push!(uplits, minimum(arr))
    @info length(arr)
end
using Plots

argmin(uplits)
plot(BLMdf.uplifts[7][1:160], dpi = 300, xlabel = "Iterations", ylabel = "Total Uplifts", label = "Bundle Level Method")
plot!(BPMdf.uplifts[7][1:160], label = "Bundle Proximal Method")
# plot!(DAdf.uplifts[7][1:160], label = "D-Adaptation")
# plot!(Polyakdf.uplifts[7], label = "True Polyak stepsize")
plot!(yscale=:log)
savefig("Stability of BPM over BLM")