using Revise
using ConvexHullPricing, DataFrames, JLD2, ProgressBars

list_of_instances = []
for file in readdir("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\ca_data"; join=true)
    push!(list_of_instances, ConvexHullPricing.Utilitaries.load_data(file))
end

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
    LP_Relax = ConvexHullPricing.Utilitaries.LP_Relaxation(instance)
    push!(BLMdf.instance, instance)
    push!(BPMdf.instance, instance)
    push!(DAdf.instance, instance)
    push!(DowGdf.instance, instance)
    push!(Polyakdf.instance, instance)
    
    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tBLM(instance, opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(BLMdf.niter, length(iter))
    push!(BLMdf.runtime, runtime)
    push!(BLMdf.oracle_time, oracle_time)
    push!(BLMdf.solve_time, solve_time)
    push!(BLMdf.relgap, relgap)
    push!(BLMdf.iterates, iter)
    push!(BLMdf.fiterates, iterf)

    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tBPM(instance, -opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(BPMdf.niter, length(iter))
    push!(BPMdf.runtime, runtime)
    push!(BPMdf.oracle_time, oracle_time)
    push!(BPMdf.solve_time, solve_time)
    push!(BPMdf.relgap, relgap)
    push!(BPMdf.iterates, iter)
    push!(BPMdf.fiterates, iterf)

    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tDA(instance, opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(DAdf.niter, length(iter))
    push!(DAdf.runtime, runtime)
    push!(DAdf.oracle_time, oracle_time)
    push!(DAdf.solve_time, solve_time)
    push!(DAdf.relgap, relgap)
    push!(DAdf.iterates, iter)
    push!(DAdf.fiterates, iterf)

    ((iter, iterf, oracle_time, solve_time), runtime) = @timed ConvexHullPricing.Optimizer.tDoWG(instance, opt, LP_Relax)
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(DowGdf.niter, length(iter))
    push!(DowGdf.runtime, runtime)
    push!(DowGdf.oracle_time, oracle_time)
    push!(DowGdf.solve_time, solve_time)
    push!(DowGdf.relgap, relgap)
    push!(DowGdf.iterates, iter)
    push!(DowGdf.fiterates, iterf)

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
optimals = load_object("results//tables//OptimalDF")
for (i, instance) in tqdm(enumerate(list_of_instances))
    opt = last(optimals[i,:].fiterates)
    benchmark(instance, opt)
end

BLMdf
BPMdf
DAdf
DowGdf
Polyakdf
save_object("results//tables//caBLMdf.jld2", BLMdf)
save_object("results//tables//caBPMdf.jld2", BPMdf)
save_object("results//tables//caDAdf.jld2", DAdf)
save_object("results//tables//caDowGdf.jld2", DowGdf)
save_object("results//tables//caPolyakdf.jld2", Polyakdf)

println("Done !!")
