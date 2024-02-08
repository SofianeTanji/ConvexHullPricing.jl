# Implementation of Dantzig-Wolfe method
using JuMP, ..Utilitaries, LinearAlgebra, Gurobi

function ColumnGeneration(instance, initial_prices, niter, eps, verbose = -1)

    # Unzip instance
    MinRunCapacity = instance.ThermalGen.MinRunCapacity
    MaxRunCapacity = instance.ThermalGen.MaxRunCapacity
    RampUp = instance.ThermalGen.RampUp
    RampDown = instance.ThermalGen.RampDown
    UpTime = instance.ThermalGen.UpTime
    DownTime = instance.ThermalGen.DownTime
    StartUp = instance.ThermalGen.StartUp
    ShutDown = instance.ThermalGen.ShutDown
    NbGen = length(MinRunCapacity)
    FixedCost = instance.ThermalGen.FixedCost
    MarginalCost = instance.ThermalGen.MarginalCost
    NoLoadConsumption = instance.ThermalGen.NoLoadConsumption
    L = instance.Load
    T = length(L)
    LostLoad = instance.LostLoad

    # Build Subproblems
    price = initial_prices
    table_subproblems = Array{Utilitaries.SubProblem}(undef, NbGen)
    for gen=1:NbGen
        model_subproblem = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
        set_silent(model_subproblem)
        set_optimizer_attribute(model_subproblem, "Threads", 8)

        # Variables
        @variable(model_subproblem, Varp[t=0:T+1], lower_bound = 0)
        @variable(model_subproblem, Varpbar[t=0:T+1], lower_bound = 0)
        @variable(model_subproblem, Varu[t=0:T+1], Bin)
        @variable(model_subproblem, Varv[t=0:T+1], Bin)
        @variable(model_subproblem, Varw[t=0:T+1], Bin)
        @variable(model_subproblem, VarCost)

        # Technical constraints
        @constraint(model_subproblem, ConstrLogical[t=1:T+1], Varu[t] - Varu[t-1] == Varv[t] - Varw[t])
        @constraint(model_subproblem, ConstrMinUpTime[t=UpTime[gen]:T+1], sum(Varv[i]  for i=t-UpTime[gen]+1:t) <= Varu[t])
        @constraint(model_subproblem, ConstrMinDownTime[t=DownTime[gen]:T+1], sum(Varw[i]  for i=t-DownTime[gen]+1:t) <= (1 - Varu[t]))
        @constraint(model_subproblem, ConstrGenLimits1[t=0:T+1], MinRunCapacity[gen]*Varu[t] <= Varp[t])
        @constraint(model_subproblem, ConstrGenLimits2[t=0:T+1], Varp[t] <= Varpbar[t])
        @constraint(model_subproblem, ConstrGenLimits3[t=0:T+1], Varpbar[t] <= MaxRunCapacity[gen]*Varu[t])
        @constraint(model_subproblem, ConstrRampUp[t=1:T+1], Varpbar[t] - Varp[t-1] <= RampUp[gen]*Varu[t-1] + StartUp[gen]*Varv[t])
        @constraint(model_subproblem, ConstrRampDown[t=1:T+1], Varpbar[t-1] - Varp[t] <= RampDown[gen]*Varu[t] + ShutDown[gen]*Varw[t])

        # Cost constraint
        @constraint(model_subproblem, VarCost - sum(NoLoadConsumption[gen] * Varu[t] + FixedCost[gen] * Varv[t] + MarginalCost[gen] * Varp[t] for t=1:T) == 0)

        # Objective
        @objective(model_subproblem, Min, sum(NoLoadConsumption[gen] * Varu[t] + FixedCost[gen] * Varv[t] + MarginalCost[gen] * Varp[t] for t=1:T) - sum(price[t]*Varp[t] for t=1:T))
        subproblem = Utilitaries.SubProblem(model_subproblem, Varp, Varpbar, Varu, Varv, Varw, VarCost, ConstrLogical, ConstrMinUpTime, ConstrMinDownTime, ConstrGenLimits1, ConstrGenLimits2, ConstrGenLimits3, ConstrRampUp, ConstrRampDown)
        table_subproblems[gen] = subproblem
    end
    if verbose > 0
        @info "[CG: Subproblems computed.]"
    end
    # Build Feasible Schedules
    ScheduleP = Dict()
    ScheduleU = Dict()
    ScheduleV = Dict()
    ScheduleW = Dict()
    ScheduleCost = Dict()
    ScheduleCounter = Int64[0 for gen=1:NbGen]

    ## Initialisation
    matching = Utilitaries.Matching(instance)
    MatchU = matching.Varu
    MatchV = matching.Varv
    MatchW = matching.Varw
    MatchP = matching.Varp

    for gen=1:NbGen
        ScheduleCounter[gen] += 1
        if NbGen > 1
            ScheduleP[gen, ScheduleCounter[gen]] = MatchP[gen, 2:T+1]
            ScheduleU[gen, ScheduleCounter[gen]] = MatchU[gen, 2:T+1]
            ScheduleV[gen, ScheduleCounter[gen]] = MatchV[gen, 2:T+1]
            ScheduleW[gen, ScheduleCounter[gen]] = MatchW[gen, 2:T+1]
            ScheduleCost[gen, ScheduleCounter[gen]] = sum(NoLoadConsumption[gen] * MatchU[gen, t] + FixedCost[gen] * MatchV[gen, t] + MarginalCost[gen] * MatchP[gen, t] for t=2:T+1)
        else
            ScheduleP[gen, ScheduleCounter[gen]] = MatchP[2:T+1]
            ScheduleU[gen, ScheduleCounter[gen]] = MatchU[2:T+1]
            ScheduleV[gen, ScheduleCounter[gen]] = MatchV[2:T+1]
            ScheduleW[gen, ScheduleCounter[gen]] = MatchW[2:T+1]
            ScheduleCost[gen, ScheduleCounter[gen]] = sum(NoLoadConsumption[gen] * MatchU[t] + FixedCost[gen] * MatchV[t] + MarginalCost[gen] * MatchP[t] for t=2:T+1)
        end
    
    end
    if verbose > 0
        @info "[CG: Initial Feasible schedules computed.]"
    end
    # Build Restricted Master Program
    restricted_model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(restricted_model)
    set_optimizer_attribute(restricted_model, "Threads", 8)
    DictZ = Dict()
    VarZ = @variable(restricted_model, [gen=1:NbGen, i=1:ScheduleCounter[gen]], lower_bound = 0, upper_bound = 1, base_name = "Z")
    for gen=1:NbGen
        for i=1:ScheduleCounter[gen]
            DictZ[gen, i] = VarZ[gen, i]
        end
    end
    VarL = @variable(restricted_model, [t=1:T], lower_bound = 0)
    @constraint(restricted_model, [t=1:T], VarL[t] <= L[t])

    ## Market Clearing Constraint
    ConstrBalance = @constraint(restricted_model, [t=1:T], sum(sum(sum(DictZ[gen,i] * prod for (time,prod) in enumerate(ScheduleP[gen,i]) if time==t) for i=1:ScheduleCounter[gen]) for gen=1:NbGen) - VarL[t] == 0)
    ConstrConvexComb = @constraint(restricted_model, [gen=1:NbGen], sum(DictZ[gen,i] for i=1:ScheduleCounter[gen]) == 1)
    # set_dual_start_value(ConstrBalance, initial_prices)
    RMP = Utilitaries.RestrictedMasterProgram(restricted_model, DictZ, VarL, ConstrBalance, ConstrConvexComb)

    ## Run the Dantzig-Wolfe iterations
    StoppingCriterion = 0
    ObjMaster = 0
    ObjVect = []
    Iterates = []
    for iter=1:niter
        if verbose > 0
            @info "[CG: Iteration $iter]"
        end
        @objective(RMP.model, Min, sum(sum(DictZ[gen, i] * ScheduleCost[gen, i] for i=1:ScheduleCounter[gen]) for gen=1:NbGen) + LostLoad * sum(L[t] - VarL[t] for t=1:T))
        optimize!(RMP.model)
        ObjMaster = objective_value(RMP.model)
        push!(ObjVect, ObjMaster)
        if iter > 1
            price = dual.(RMP.ConstrBalance)
            push!(Iterates, price)
        else
            price = initial_prices
            push!(Iterates, price)
        end
        PiDual = dual.(RMP.ConstrConvexComb)
        StoppingCriterion = 1
        for gen=1:NbGen
            SubProblem = table_subproblems[gen]
            @objective(SubProblem.model, Min, sum(NoLoadConsumption[gen] * SubProblem.Varu[t] + FixedCost[gen] * SubProblem.Varv[t] + MarginalCost[gen] * SubProblem.Varp[t] for t=1:T) - sum(price[t] * SubProblem.Varp[t] for t=1:T))
            optimize!(SubProblem.model)
            termination_status(SubProblem.model)
            ReducedCost = objective_value(SubProblem.model) - PiDual[gen]
            if ReducedCost <= -eps
                StoppingCriterion = 0
                SubP = value.(SubProblem.Varp).data
                SubU = value.(SubProblem.Varu).data
                SubV = value.(SubProblem.Varv).data
                SubW = value.(SubProblem.Varw).data

                ScheduleCounter[gen] += 1
                ScheduleP[gen, ScheduleCounter[gen]] = SubP[2:T+1]
                ScheduleU[gen, ScheduleCounter[gen]] = SubP[2:T+1]
                ScheduleV[gen, ScheduleCounter[gen]] = SubP[2:T+1]
                ScheduleW[gen, ScheduleCounter[gen]] = SubP[2:T+1]
                ScheduleCost[gen, ScheduleCounter[gen]] = sum(NoLoadConsumption[gen] * SubU[t] + FixedCost[gen] * SubV[t] + MarginalCost[gen] * SubP[t] for t=2:T+1)

                NewVarZ = @variable(RMP.model, [[gen], [ScheduleCounter[gen]]], lower_bound = 0, upper_bound = 1, base_name = "Z")
                RMP.VarZ[gen, ScheduleCounter[gen]] = NewVarZ[gen, ScheduleCounter[gen]]
                
                for t=1:T
                    set_normalized_coefficient(RMP.ConstrBalance[t], RMP.VarZ[gen, ScheduleCounter[gen]], SubP[t+1])
                end
                set_normalized_coefficient(RMP.ConstrConvexComb[gen], RMP.VarZ[gen, ScheduleCounter[gen]], 1)
            end
        end
        if StoppingCriterion == 1
            break
        end
    end
    f_iterates = Float64[]
    for ρ in Iterates
        push!(f_iterates, Utilitaries.exact_oracle(instance, ρ)[1])
    end
    return last(Iterates), Iterates, f_iterates
end

function tColumnGeneration(instance, initial_prices, budget, eps, verbose = -1)

    # Unzip instance
    MinRunCapacity = instance.ThermalGen.MinRunCapacity
    MaxRunCapacity = instance.ThermalGen.MaxRunCapacity
    RampUp = instance.ThermalGen.RampUp
    RampDown = instance.ThermalGen.RampDown
    UpTime = instance.ThermalGen.UpTime
    DownTime = instance.ThermalGen.DownTime
    StartUp = instance.ThermalGen.StartUp
    ShutDown = instance.ThermalGen.ShutDown
    NbGen = length(MinRunCapacity)
    FixedCost = instance.ThermalGen.FixedCost
    MarginalCost = instance.ThermalGen.MarginalCost
    NoLoadConsumption = instance.ThermalGen.NoLoadConsumption
    L = instance.Load
    T = length(L)
    LostLoad = instance.LostLoad
    it_time = @elapsed begin
    # Build Subproblems
    price = initial_prices
    table_subproblems = Array{Utilitaries.SubProblem}(undef, NbGen)
    for gen=1:NbGen
        model_subproblem = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
        set_silent(model_subproblem)
        set_optimizer_attribute(model_subproblem, "Threads", 8)

        # Variables
        @variable(model_subproblem, Varp[t=0:T+1], lower_bound = 0)
        @variable(model_subproblem, Varpbar[t=0:T+1], lower_bound = 0)
        @variable(model_subproblem, Varu[t=0:T+1], Bin)
        @variable(model_subproblem, Varv[t=0:T+1], Bin)
        @variable(model_subproblem, Varw[t=0:T+1], Bin)
        @variable(model_subproblem, VarCost)

        # Technical constraints
        @constraint(model_subproblem, ConstrLogical[t=1:T+1], Varu[t] - Varu[t-1] == Varv[t] - Varw[t])
        @constraint(model_subproblem, ConstrMinUpTime[t=UpTime[gen]:T+1], sum(Varv[i]  for i=t-UpTime[gen]+1:t) <= Varu[t])
        @constraint(model_subproblem, ConstrMinDownTime[t=DownTime[gen]:T+1], sum(Varw[i]  for i=t-DownTime[gen]+1:t) <= (1 - Varu[t]))
        @constraint(model_subproblem, ConstrGenLimits1[t=0:T+1], MinRunCapacity[gen]*Varu[t] <= Varp[t])
        @constraint(model_subproblem, ConstrGenLimits2[t=0:T+1], Varp[t] <= Varpbar[t])
        @constraint(model_subproblem, ConstrGenLimits3[t=0:T+1], Varpbar[t] <= MaxRunCapacity[gen]*Varu[t])
        @constraint(model_subproblem, ConstrRampUp[t=1:T+1], Varpbar[t] - Varp[t-1] <= RampUp[gen]*Varu[t-1] + StartUp[gen]*Varv[t])
        @constraint(model_subproblem, ConstrRampDown[t=1:T+1], Varpbar[t-1] - Varp[t] <= RampDown[gen]*Varu[t] + ShutDown[gen]*Varw[t])

        # Cost constraint
        @constraint(model_subproblem, VarCost - sum(NoLoadConsumption[gen] * Varu[t] + FixedCost[gen] * Varv[t] + MarginalCost[gen] * Varp[t] for t=1:T) == 0)

        # Objective
        @objective(model_subproblem, Min, sum(NoLoadConsumption[gen] * Varu[t] + FixedCost[gen] * Varv[t] + MarginalCost[gen] * Varp[t] for t=1:T) - sum(price[t]*Varp[t] for t=1:T))
        subproblem = Utilitaries.SubProblem(model_subproblem, Varp, Varpbar, Varu, Varv, Varw, VarCost, ConstrLogical, ConstrMinUpTime, ConstrMinDownTime, ConstrGenLimits1, ConstrGenLimits2, ConstrGenLimits3, ConstrRampUp, ConstrRampDown)
        table_subproblems[gen] = subproblem
    end
    if verbose > 0
        @info "[CG: Subproblems computed.]"
    end
    # Build Feasible Schedules
    ScheduleP = Dict()
    ScheduleU = Dict()
    ScheduleV = Dict()
    ScheduleW = Dict()
    ScheduleCost = Dict()
    ScheduleCounter = Int64[0 for gen=1:NbGen]
    if verbose > 0
        @info "[CG: Computing initial feasible schedules.]"
    end
    ## Initialisation
    matching = Utilitaries.Matching(instance)
    MatchU = matching.Varu
    MatchV = matching.Varv
    MatchW = matching.Varw
    MatchP = matching.Varp

    for gen=1:NbGen
        if verbose > 0
            @info "[CG: Feasible schedule initial computation, generator $gen]"
        end
        ScheduleCounter[gen] += 1
        if NbGen > 1
            ScheduleP[gen, ScheduleCounter[gen]] = MatchP[gen, 2:T+1]
            ScheduleU[gen, ScheduleCounter[gen]] = MatchU[gen, 2:T+1]
            ScheduleV[gen, ScheduleCounter[gen]] = MatchV[gen, 2:T+1]
            ScheduleW[gen, ScheduleCounter[gen]] = MatchW[gen, 2:T+1]
            ScheduleCost[gen, ScheduleCounter[gen]] = sum(NoLoadConsumption[gen] * MatchU[gen, t] + FixedCost[gen] * MatchV[gen, t] + MarginalCost[gen] * MatchP[gen, t] for t=2:T+1)
        else
            ScheduleP[gen, ScheduleCounter[gen]] = MatchP[2:T+1]
            ScheduleU[gen, ScheduleCounter[gen]] = MatchU[2:T+1]
            ScheduleV[gen, ScheduleCounter[gen]] = MatchV[2:T+1]
            ScheduleW[gen, ScheduleCounter[gen]] = MatchW[2:T+1]
            ScheduleCost[gen, ScheduleCounter[gen]] = sum(NoLoadConsumption[gen] * MatchU[t] + FixedCost[gen] * MatchV[t] + MarginalCost[gen] * MatchP[t] for t=2:T+1)
        end
    
    end
    if verbose > 0
        @info "[CG: Feasible schedules computed.]"
    end
    # Build Restricted Master Program
    restricted_model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(restricted_model)
    set_optimizer_attribute(restricted_model, "Threads", 8)
    DictZ = Dict()
    VarZ = @variable(restricted_model, [gen=1:NbGen, i=1:ScheduleCounter[gen]], lower_bound = 0, upper_bound = 1, base_name = "Z")
    for gen=1:NbGen
        for i=1:ScheduleCounter[gen]
            DictZ[gen, i] = VarZ[gen, i]
        end
    end
    VarL = @variable(restricted_model, [t=1:T], lower_bound = 0)
    @constraint(restricted_model, [t=1:T], VarL[t] <= L[t])

    ## Market Clearing Constraint
    ConstrBalance = @constraint(restricted_model, [t=1:T], sum(sum(sum(DictZ[gen,i] * prod for (time,prod) in enumerate(ScheduleP[gen,i]) if time==t) for i=1:ScheduleCounter[gen]) for gen=1:NbGen) - VarL[t] == 0)
    ConstrConvexComb = @constraint(restricted_model, [gen=1:NbGen], sum(DictZ[gen,i] for i=1:ScheduleCounter[gen]) == 1)
    # set_dual_start_value(ConstrBalance, initial_prices)
    RMP = Utilitaries.RestrictedMasterProgram(restricted_model, DictZ, VarL, ConstrBalance, ConstrConvexComb)

    ## Run the Dantzig-Wolfe iterations
    StoppingCriterion = 0
    ObjMaster = 0
    ObjVect = []
    Iterates = []
    end
    time_vector = [it_time]
    iter = 1
    while time_vector[end] <= budget
        if verbose > 0
            @info "[CG: Iteration $iter]"
        end
        it_time = @elapsed begin
        @objective(RMP.model, Min, sum(sum(DictZ[gen, i] * ScheduleCost[gen, i] for i=1:ScheduleCounter[gen]) for gen=1:NbGen) + LostLoad * sum(L[t] - VarL[t] for t=1:T))
        optimize!(RMP.model)
        ObjMaster = objective_value(RMP.model)
        push!(ObjVect, ObjMaster)
        if iter > 1
            price = dual.(RMP.ConstrBalance)
            push!(Iterates, price)
        else
            price = initial_prices
            push!(Iterates, price)
        end
        PiDual = dual.(RMP.ConstrConvexComb)
        StoppingCriterion = 1
        for gen=1:NbGen
            SubProblem = table_subproblems[gen]
            @objective(SubProblem.model, Min, sum(NoLoadConsumption[gen] * SubProblem.Varu[t] + FixedCost[gen] * SubProblem.Varv[t] + MarginalCost[gen] * SubProblem.Varp[t] for t=1:T) - sum(price[t] * SubProblem.Varp[t] for t=1:T))
            optimize!(SubProblem.model)
            termination_status(SubProblem.model)
            ReducedCost = objective_value(SubProblem.model) - PiDual[gen]
            if ReducedCost <= eps
                StoppingCriterion = 0
                SubP = value.(SubProblem.Varp).data
                SubU = value.(SubProblem.Varu).data
                SubV = value.(SubProblem.Varv).data
                SubW = value.(SubProblem.Varw).data

                ScheduleCounter[gen] += 1
                ScheduleP[gen, ScheduleCounter[gen]] = SubP[2:T+1]
                ScheduleU[gen, ScheduleCounter[gen]] = SubP[2:T+1]
                ScheduleV[gen, ScheduleCounter[gen]] = SubP[2:T+1]
                ScheduleW[gen, ScheduleCounter[gen]] = SubP[2:T+1]
                ScheduleCost[gen, ScheduleCounter[gen]] = sum(NoLoadConsumption[gen] * SubU[t] + FixedCost[gen] * SubV[t] + MarginalCost[gen] * SubP[t] for t=2:T+1)

                NewVarZ = @variable(RMP.model, [[gen], [ScheduleCounter[gen]]], lower_bound = 0, upper_bound = 1, base_name = "Z")
                RMP.VarZ[gen, ScheduleCounter[gen]] = NewVarZ[gen, ScheduleCounter[gen]]
                
                for t=1:T
                    set_normalized_coefficient(RMP.ConstrBalance[t], RMP.VarZ[gen, ScheduleCounter[gen]], SubP[t+1])
                end
                set_normalized_coefficient(RMP.ConstrConvexComb[gen], RMP.VarZ[gen, ScheduleCounter[gen]], 1)
            end
        end
        if StoppingCriterion == 1
            break
        end
        end
        push!(time_vector, it_time + time_vector[end])
        iter += 1
    end
    f_iterates = Float64[]
    for ρ in Iterates
        push!(f_iterates, Utilitaries.exact_oracle(instance, ρ)[1])
    end
    return last(Iterates), Iterates, f_iterates, time_vector
end