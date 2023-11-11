using Revise
using ConvexHullPricing, DataFrames, JLD2, ProgressBars

const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer

@info "Loading instances ..."
list_of_instances = []
for file in readdir("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\data\\belgian"; join=true)
    push!(list_of_instances, UT.load_data(file))
end
@info "BE instances are loaded !"

@info "Loading instances ..."
for file in readdir("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\data\\ca"; join=true)
    # push!(list_of_instances, UT.load_data(file))
end
@info "CA instances are loaded !"
list_of_instances
CGdf = DataFrame(
        instance = [],
        niter = [],
        runtime = [],
        oracle_time = [],
        solve_time = [],
        relgap = [],
        iterates = [],
        fiterates = []
)
@info "Dataframe for CG is ready to be filled !"

function benchmark(instance, opt)
    @info "Computing initial iterate"
    LP_Relax = UT.LP_Relaxation(instance)
    push!(CGdf.instance, instance)
    
    @info "Running the CG method"
    ((iter, solve_time), runtime) = @timed OPT.tCG(instance, LP_Relax)
    iterf = []
    for x in iter
        push!(iterf, UT.fast_oracle(instance, x)[1])
    end
    relgap = abs(abs(last(iterf)) - abs(opt))/maximum([abs(last(iterf)), abs(opt)])
    push!(CGdf.niter, length(iter))
    push!(CGdf.runtime, runtime)
    push!(CGdf.oracle_time, 0)
    push!(CGdf.solve_time, solve_time)
    push!(CGdf.relgap, relgap)
    push!(CGdf.iterates, iter)
    push!(CGdf.fiterates, iterf)
end

beoptimals = load_object("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\results\\objects\\beOptimalDF.jld2")
# caoptimals = load_object("C:\\Users\\Sofiane\\Desktop\\ConvexHullPricing\\results\\objects\\caOptimalDF.jld2")
# optimals = vcat(beoptimals, caoptimals)
@info "Optimal values of each instance loaded."
@info "Running benchmarks"

@info "Max iter = 500."
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
    Load = [6 11 16 11],
    ThermalGen = ThermalGen
)
LP_Relax = UT.LP_Relaxation(instance)
@timed OPT.tCG(instance, LP_Relax)

for (i, instance) in tqdm(enumerate(list_of_instances))
    opt = last(beoptimals[i,:].fiterates)
    benchmark(instance, opt)
end
@info "Benchmarks done ! Now showing df"
@show CGdf
println("Everything is done !!")
