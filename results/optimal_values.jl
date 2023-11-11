using ConvexHullPricing, JLD2, ProgressBars, DataFrames
df = DataFrame(
    instance = [],
    relgap = [],
    iterates = [],
    fiterates = []
)
@info "Loading instances ..."
list_of_instances = []
for file in readdir("\\data\\ca"; join=true)
    push!(list_of_instances, ConvexHullPricing.Utilitaries.load_data(file))
end
@info "CA instances loaded !"
for (i, instance) in tqdm(enumerate(list_of_instances))
    push!(df.instance, instance)
    LP_Relax = ConvexHullPricing.Utilitaries.LP_Relaxation(instance)
    xstar, iter, fiter, relgap = ConvexHullPricing.Optimizer.GetOptimal(instance, LP_Relax, 1200, 1e-7, .98)
    @info "For instance #$(i), relative gap reached = $(relgap)"
    push!(df.relgap, relgap)
    push!(df.iterates, iter)
    push!(df.fiterates, fiter)
end

for i=1:20
    @info length(df.iterates[i])
end

save_object("results//tables//caOptimalDF.jld2", df)


# Now for Belgian

bedf = DataFrame(
    instance = [],
    relgap = [],
    iterates = [],
    fiterates = []
)
@info "Here we go !"
list_of_instances = []
for file in readdir("\\data\\belgian"; join=true)
    push!(list_of_instances, ConvexHullPricing.Utilitaries.load_data(file))
end
@info "Here we go !"
for (i, instance) in tqdm(enumerate(list_of_instances))
    push!(bedf.instance, instance)
    LP_Relax = ConvexHullPricing.Utilitaries.LP_Relaxation(instance)
    xstar, iter, fiter, relgap = ConvexHullPricing.Optimizer.GetOptimal(instance, LP_Relax, 1200, 1e-7, .98)
    @info "For instance #$(i), relative gap reached = $(relgap)"
    push!(bedf.relgap, relgap)
    push!(bedf.iterates, iter)
    push!(bedf.fiterates, fiter)
end

for i=1:8
    @info (bedf.relgap[i])
    @info (length(bedf.iterates[i]))
end

save_object("results//tables//beOptimalDF.jld2", bedf)