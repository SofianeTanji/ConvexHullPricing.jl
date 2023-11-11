using ConvexHullPricing
using Gurobi
using Plots

println("Here we go !")
list_of_instances = []
for file in readdir("\\data\\ca"; join=true)
    push!(list_of_instances, ConvexHullPricing.Utilitaries.load_data(file))
end
println("Instances loaded")

list_stepsizes = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]

distances = zeros(length(list_of_instances), length(list_stepsizes))

df = load_object("results//tables//CaDF.jld2")

for (i, instance) in tqdm(enumerate(list_of_instances))
    LP_Relax = ConvexHullPricing.Utilitaries.LP_Relaxation(instance)
    opt = last(df[i,:].fiterates)
    for (j, step) in tqdm(enumerate(list_stepsizes))
        xstar, iter, fiter = ConvexHullPricing.Optimizer.BundleProximalMethod(instance, LP_Relax, 5, .98, step)
        distances[i, j] = abs(last(fiter) - opt)/maximum([opt, last(fiter)])
    end
    println("For instance #$(i), best relative gap reached = $(minimum(distances[i,:]))")
end

for i=1:20
    @show maximum(distances[i,:])
end