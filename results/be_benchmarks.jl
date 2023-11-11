using Revise
using ConvexHullPricing, DataFrames, JLD2, ProgressBars

@info "Loading instances ..."
list_of_instances = []
for file in readdir("\\data\\belgian"; join=true)
    push!(list_of_instances, ConvexHullPricing.Utilitaries.load_data(file))
end
@info "BE instances are loaded !"

function init()
    BLMdf = DataFrame(
        instance = [],
        niter = [],
        runtime = [],
        oracle_time = [],
        solve_time = [],
        relgap = [],
        iterates = [],
        fiterates = []
    )
    BPMdf = DataFrame(
        instance = [],
        niter = [],
        runtime = [],
        oracle_time = [],
        solve_time = [],
        relgap = [],
        iterates = [],
        fiterates = []
    )

    DAdf = DataFrame(
        instance = [],
        niter = [],
        runtime = [],
        oracle_time = [],
        solve_time = [],
        relgap = [],
        iterates = [],
        fiterates = []
    )

    DowGdf = DataFrame(
        instance = [],
        niter = [],
        runtime = [],
        oracle_time = [],
        solve_time = [],
        relgap = [],
        iterates = [],
        fiterates = []
    )

    Polyakdf = DataFrame(
        instance = [],
        niter = [],
        runtime = [],
        oracle_time = [],
        solve_time = [],
        relgap = [],
        iterates = [],
        fiterates = []
    )
    return BLMdf, BPMdf, DAdf, DowGdf, Polyakdf
end
function benchmark(instance, opt)
    @info "Computing initial iterate"
    LP_Relax = ConvexHullPricing.Utilitaries.LP_Relaxation(instance)
    push!(BLMdf.instance, instance)
    push!(BPMdf.instance, instance)
    push!(DAdf.instance, instance)
    push!(DowGdf.instance, instance)
    push!(Polyakdf.instance, instance)
    
    @info "Running the BLM method"
    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tBLM(instance, opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(BLMdf.niter, length(iter))
    push!(BLMdf.runtime, runtime)
    push!(BLMdf.oracle_time, oracle_time)
    push!(BLMdf.solve_time, solve_time)
    push!(BLMdf.relgap, relgap)
    push!(BLMdf.iterates, iter)
    push!(BLMdf.fiterates, iterf)

    @info "Running the BPM method"
    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tBPM(instance, -opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(BPMdf.niter, length(iter))
    push!(BPMdf.runtime, runtime)
    push!(BPMdf.oracle_time, oracle_time)
    push!(BPMdf.solve_time, solve_time)
    push!(BPMdf.relgap, relgap)
    push!(BPMdf.iterates, iter)
    push!(BPMdf.fiterates, iterf)

    @info "Running the DAdaptation method"
    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tDA(instance, opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(DAdf.niter, length(iter))
    push!(DAdf.runtime, runtime)
    push!(DAdf.oracle_time, oracle_time)
    push!(DAdf.solve_time, solve_time)
    push!(DAdf.relgap, relgap)
    push!(DAdf.iterates, iter)
    push!(DAdf.fiterates, iterf)

    @info "Running the DowG method"
    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tDoWG(instance, opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(DowGdf.niter, length(iter))
    push!(DowGdf.runtime, runtime)
    push!(DowGdf.oracle_time, oracle_time)
    push!(DowGdf.solve_time, solve_time)
    push!(DowGdf.relgap, relgap)
    push!(DowGdf.iterates, iter)
    push!(DowGdf.fiterates, iterf)

    @info "Running the Polyak method"
    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tPolyak(instance, opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(Polyakdf.niter, length(iter))
    push!(Polyakdf.runtime, runtime)
    push!(Polyakdf.oracle_time, oracle_time)
    push!(Polyakdf.solve_time, solve_time)
    push!(Polyakdf.relgap, relgap)
    push!(Polyakdf.iterates, iter)
    push!(Polyakdf.fiterates, iterf)
end

BLMdf, BPMdf, DAdf, DowGdf, Polyakdf = init()
@info "Dataframes for each method are ready to be filled !"

optimals = load_object("results//tables//beOptimalDF.jld2")

@info "Optimal values of each instance loaded."

@info "Running benchmarks for each method."
@info "Max iter = 500. Good hyperparameters guessed but not optimized."
for (i, instance) in tqdm(enumerate(list_of_instances))
    opt = last(optimals[i,:].fiterates)
    benchmark(instance, opt)
end
@info "Benchmarks done !"

@show BLMdf
@show BPMdf
@show DAdf
@show DowGdf
@show Polyakdf

save_object("results//tables//beBLMdf.jld2", BLMdf)
save_object("results//tables//beBPMdf.jld2", BPMdf)
save_object("results//tables//beDAdf.jld2", DAdf)
save_object("results//tables//beDowGdf.jld2", DowGdf)
save_object("results//tables//bePolyakdf.jld2", Polyakdf)

println("Everything is done !!")

# ---------- # Correct Relative Gaps
BLMdf, BPMdf, DAdf, DowGdf, Polyakdf = load_object("results//tables//beBLMdf.jld2"), load_object("results//tables//beBPMdf.jld2"), load_object("results//tables//beDAdf.jld2"), load_object("results//tables//beDowGdf.jld2"), load_object("results//tables//bePolyakdf.jld2")
optimals = load_object("results//tables//beOptimalDF.jld2")
for i=1:8
    opt = optimals.fiterates[i][end]
    if length(DAdf.fiterates[i]) > 1
        BLMdf.relgap[i] = abs(BLMdf.fiterates[i][end - 1] - optimals.fiterates[i][end])/abs(optimals.fiterates[i][end])
    end
end

@show BLMdf
@show BPMdf
@show DAdf
@show DowGdf
@show Polyakdf

save_object("results//tables//beBLMdf.jld2", BLMdf)
save_object("results//tables//beBPMdf.jld2", BPMdf)
save_object("results//tables//beDAdf.jld2", DAdf)
save_object("results//tables//beDowGdf.jld2", DowGdf)
save_object("results//tables//bePolyakdf.jld2", Polyakdf)

@info "Finito les benchmarks"