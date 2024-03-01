# Implementation of defined oracles in (Tanji et al., 2023)
using JuMP, Gurobi, LinearAlgebra, Random

function Matching(instance)
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

    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model)
    set_optimizer_attributes(model, "MIPGap" => 1e-3)

    @variable(model, Varp[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varpbar[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varu[g = 1:NbGen, t = 0:T+1], Bin)
    @variable(model, Varv[g = 1:NbGen, t = 0:T+1], Bin)
    @variable(model, Varw[g = 1:NbGen, t = 0:T+1], Bin)

    @variable(model, VarL[t = 1:T], lower_bound = 0)
    @constraint(model, [t = 1:T], VarL[t] <= L[t])

    @constraint(model, loads[t = 1:T], sum(Varp[g, t] for g = 1:NbGen) == VarL[t]) # Balance constraint

    @constraint(
        model,
        ConstrLogical[gen = 1:NbGen, t = 1:T+1],
        Varu[gen, t] - Varu[gen, t-1] == Varv[gen, t] - Varw[gen, t]
    )
    @constraint(
        model,
        ConstrMinUpTime[gen = 1:NbGen, t = UpTime[gen]:T+1],
        sum(Varv[gen, i] for i = t-UpTime[gen]+1:t) <= Varu[gen, t]
    )
    @constraint(
        model,
        ConstrMinDownTime[gen = 1:NbGen, t = DownTime[gen]:T+1],
        sum(Varw[gen, i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[gen, t])
    )
    @constraint(
        model,
        ConstrGenLimits1[gen = 1:NbGen, t = 0:T+1],
        MinRunCapacity[gen] * Varu[gen, t] <= Varp[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits2[gen = 1:NbGen, t = 0:T+1],
        Varp[gen, t] <= Varpbar[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits3[gen = 1:NbGen, t = 0:T+1],
        Varpbar[gen, t] <= MaxRunCapacity[gen] * Varu[gen, t]
    )
    @constraint(
        model,
        ConstrRampUp[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t] - Varp[gen, t-1] <=
        RampUp[gen] * Varu[gen, t-1] + StartUp[gen] * Varv[gen, t]
    )
    @constraint(
        model,
        ConstrRampDown[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t-1] - Varp[gen, t] <=
        RampDown[gen] * Varu[gen, t] + ShutDown[gen] * Varw[gen, t]
    )

    @objective(
        model,
        Min,
        sum(
            sum(
                NoLoadConsumption[gen] * Varu[gen, t] +
                FixedCost[gen] * Varv[gen, t] +
                MarginalCost[gen] * Varp[gen, t] for t = 1:T
            ) for gen = 1:NbGen
        ) + LostLoad * sum(L[t] - VarL[t] for t = 1:T)
    )
    optimize!(model)

    ValU, ValV, ValW, ValP, ValPbar, ValL = value.(Varu).data,
    value.(Varv).data,
    value.(Varw).data,
    value.(Varp).data,
    value.(Varpbar).data,
    value.(VarL)
    return MatchingSolution(objective_value(model), ValU, ValV, ValW, ValP, ValPbar, ValL)
end

function GetMShift(instance)
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

    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model)

    @variable(model, Varp[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varpbar[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varu[g = 1:NbGen, t = 0:T+1], Bin)
    @variable(model, Varv[g = 1:NbGen, t = 0:T+1], Bin)
    @variable(model, Varw[g = 1:NbGen, t = 0:T+1], Bin)

    @variable(model, VarL[t = 1:T], lower_bound = 0)
    @constraint(model, [t = 1:T], VarL[t] <= L[t])

    @constraint(model, loads[t = 1:T], sum(Varp[g, t] for g = 1:NbGen) == VarL[t]) # Balance constraint

    @constraint(
        model,
        ConstrLogical[gen = 1:NbGen, t = 1:T+1],
        Varu[gen, t] - Varu[gen, t-1] == Varv[gen, t] - Varw[gen, t]
    )
    @constraint(
        model,
        ConstrMinUpTime[gen = 1:NbGen, t = UpTime[gen]:T+1],
        sum(Varv[gen, i] for i = t-UpTime[gen]+1:t) <= Varu[gen, t]
    )
    @constraint(
        model,
        ConstrMinDownTime[gen = 1:NbGen, t = DownTime[gen]:T+1],
        sum(Varw[gen, i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[gen, t])
    )
    @constraint(
        model,
        ConstrGenLimits1[gen = 1:NbGen, t = 0:T+1],
        MinRunCapacity[gen] * Varu[gen, t] <= Varp[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits2[gen = 1:NbGen, t = 0:T+1],
        Varp[gen, t] <= Varpbar[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits3[gen = 1:NbGen, t = 0:T+1],
        Varpbar[gen, t] <= MaxRunCapacity[gen] * Varu[gen, t]
    )
    @constraint(
        model,
        ConstrRampUp[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t] - Varp[gen, t-1] <=
        RampUp[gen] * Varu[gen, t-1] + StartUp[gen] * Varv[gen, t]
    )
    @constraint(
        model,
        ConstrRampDown[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t-1] - Varp[gen, t] <=
        RampDown[gen] * Varu[gen, t] + ShutDown[gen] * Varw[gen, t]
    )

    @objective(
        model,
        Min,
        sum(
            sum(
                NoLoadConsumption[gen] * Varu[gen, t] +
                FixedCost[gen] * Varv[gen, t] +
                MarginalCost[gen] * Varp[gen, t] for t = 1:T
            ) for gen = 1:NbGen
        ) + LostLoad * sum(L[t] - VarL[t] for t = 1:T)
    )
    optimize!(model)

    U = value.(model[:Varu]).data
    V = value.(model[:Varv]).data
    W = value.(model[:Varw]).data
    P = value.(model[:Varp]).data
    Pbar = value.(model[:Varpbar]).data
    LL = value.(model[:VarL])
    return sum(P[t]^2 + Pbar[t]^2 + U[t]^2 + V[t]^2 + W[t]^2 + LL[t]^2 for t = 1:T)
end

function set_T(UT, T_max, nb_gen)
    index_gen = zeros(nb_gen)
    for g = 1:nb_gen
        index = 0
        # Intervals for generator "on" on prior, intervals [0,b]
        # index+=T_max+1; # [0,b] for b=0:T_max
        # Intervals [a,b] such that 1 <= a <= a+(UT-1) <= b <= T_max
        for a = 1:T_max
            for b = a+UT[g]-1:T_max
                if (1 <= a && a <= a + UT[g] - 1 && a + UT[g] - 1 <= b && b <= T_max)
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
    A = zeros(nb_gen, convert(Int64, maximum(index_gen)))
    B = zeros(nb_gen, convert(Int64, maximum(index_gen)))

    for g = 1:nb_gen
        index = 1
        for a = 1:T_max
            for b = a+UT[g]-1:T_max
                if (1 <= a && a <= a + UT[g] - 1 && a + UT[g] - 1 <= b && b <= T_max)
                    A[g, index] = a
                    B[g, index] = b
                    index += 1
                end
            end
        end
    end
    A = convert(Matrix{Int64}, A)
    B = convert(Matrix{Int64}, B)
    index_gen = convert(Array{Int64}, index_gen)
    return (A, B, index_gen)
end

function ExtendedFormulation(instance)
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

    (A, B, NbIntervalsGen) = set_T(UpTime, T, NbGen)

    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model)

    # Variablesg=1:nb_gen, i=1:nb_intervals_gen[g], t=0:T_max+1
    @variable(model, Varp[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = 0:T+1], lower_bound = 0)
    @variable(
        model,
        Varpbar[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = 0:T+1],
        lower_bound = 0
    )
    @variable(model, VarpTime[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, VarpbarTime[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Vargamma[g = 1:NbGen, i = 1:NbIntervalsGen[g]], lower_bound = 0)
    @variable(model, Varu[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)
    @variable(model, Varv[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)
    @variable(model, Varw[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)
    @variable(model, VarL[t = 1:T], lower_bound = 0)
    @constraint(model, [t = 1:T], VarL[t] <= L[t])

    @objective(
        model,
        Min,
        sum(
            sum(
                NoLoadConsumption[gen] * Varu[gen, t] +
                FixedCost[gen] * Varv[gen, t] +
                MarginalCost[gen] * VarpTime[gen, t] for t = 1:T
            ) for gen = 1:NbGen
        ) + LostLoad * sum(L[t] - VarL[t] for t = 1:T)
    )

    @constraint(model, loads[t = 1:T], sum(VarpTime[g, t] for g = 1:NbGen) == VarL[t]) # Balance constraint

    # Feasible Dispatch Polytope
    @constraint(
        model,
        MinOutput[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]:B[g, i]],
        -Varp[g, i, t] <= -Vargamma[g, i] * MinRunCapacity[g]
    )
    @constraint(
        model,
        MaxOutput[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]:B[g, i]],
        Varp[g, i, t] - Varpbar[g, i, t] <= 0
    )

    @constraint(
        model,
        UpperBoundOnPbar1[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]:B[g, i]],
        Varpbar[g, i, t] <= Vargamma[g, i] * MaxRunCapacity[g]
    )
    @constraint(
        model,
        UpperBoundOnPbar2[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]:B[g, i]],
        Varpbar[g, i, t] <= Vargamma[g, i] * (StartUp[g] + (t - A[g, i]) * RampUp[g])
    )
    @constraint(
        model,
        UpperBoundOnPbar3[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]:B[g, i]],
        Varpbar[g, i, t] <= Vargamma[g, i] * (ShutDown[g] + (B[g, i] - t) * RampDown[g])
    )

    @constraint(
        model,
        LimitPowerJumpsUp1[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]+1:B[g, i]],
        Varpbar[g, i, t] - Varp[g, i, t-1] <= Vargamma[g, i] * RampUp[g]
    )
    @constraint(
        model,
        LimitPowerJumpsUp2[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]+1:B[g, i]],
        Varpbar[g, i, t] - Varp[g, i, t-1] <=
        Vargamma[g, i] * (ShutDown[g] + (B[g, i] - t) * RampDown[g] - MinRunCapacity[g])
    )

    @constraint(
        model,
        LimitPowerJumpsDown1[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]+1:B[g, i]],
        Varpbar[g, i, t-1] - Varp[g, i, t] <= Vargamma[g, i] * RampDown[g]
    )
    @constraint(
        model,
        LimitPowerJumpsDown2[g = 1:NbGen, i = 1:NbIntervalsGen[g], t = A[g, i]+1:B[g, i]],
        Varpbar[g, i, t-1] - Varp[g, i, t] <=
        Vargamma[g, i] * (StartUp[g] + (t - A[g, i]) * RampUp[g] - MinRunCapacity[g])
    )

    # Packing Dispatch Polytope
    @constraint(
        model,
        MinDownTime[g = 1:NbGen, t = 1:T],
        sum(
            Vargamma[g, i] for
            i = 1:NbIntervalsGen[g] if (A[g, i] <= t && t <= B[g, i] + DownTime[g])
        ) <= 1
    )
    @constraint(
        model,
        ProdTime[g = 1:NbGen, t = 1:T],
        VarpTime[g, t] == sum(Varp[g, i, t] for i = 1:NbIntervalsGen[g])
    )
    @constraint(
        model,
        ProdBarTime[g = 1:NbGen, t = 1:T],
        VarpbarTime[g, t] == sum(Varp[g, i, t] for i = 1:NbIntervalsGen[g])
    )

    @constraint(
        model,
        CommitStatus[g = 1:NbGen, t = 0:T+1],
        sum(Vargamma[g, i] for i = 1:NbIntervalsGen[g] if (A[g, i] <= t && t <= B[g, i])) ==
        Varu[g, t]
    )
    @constraint(
        model,
        StartUpStatus[g = 1:NbGen, t = 0:T+1],
        sum(Vargamma[g, i] for i = 1:NbIntervalsGen[g] if (A[g, i] == t)) == Varv[g, t]
    )
    @constraint(
        model,
        ShutDownStatus[g = 1:NbGen, t = 0:T+1],
        sum(Vargamma[g, i] for i = 1:NbIntervalsGen[g] if (B[g, i] + 1 == t)) == Varw[g, t]
    )

    # Optimize !
    optimize!(model)

    return shadow_price.(loads)
end

function LP_Relaxation(instance)
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

    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model)

    @variable(model, Varp[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varpbar[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varu[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)
    @variable(model, Varv[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)
    @variable(model, Varw[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)

    @variable(model, VarL[t = 1:T], lower_bound = 0)
    @constraint(model, [t = 1:T], VarL[t] <= L[t])

    @constraint(model, loads[t = 1:T], sum(Varp[g, t] for g = 1:NbGen) == VarL[t]) # Balance constraint

    @constraint(
        model,
        ConstrLogical[gen = 1:NbGen, t = 1:T+1],
        Varu[gen, t] - Varu[gen, t-1] == Varv[gen, t] - Varw[gen, t]
    )
    @constraint(
        model,
        ConstrMinUpTime[gen = 1:NbGen, t = UpTime[gen]:T+1],
        sum(Varv[gen, i] for i = t-UpTime[gen]+1:t) <= Varu[gen, t]
    )
    @constraint(
        model,
        ConstrMinDownTime[gen = 1:NbGen, t = DownTime[gen]:T+1],
        sum(Varw[gen, i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[gen, t])
    )
    @constraint(
        model,
        ConstrGenLimits1[gen = 1:NbGen, t = 0:T+1],
        MinRunCapacity[gen] * Varu[gen, t] <= Varp[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits2[gen = 1:NbGen, t = 0:T+1],
        Varp[gen, t] <= Varpbar[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits3[gen = 1:NbGen, t = 0:T+1],
        Varpbar[gen, t] <= MaxRunCapacity[gen] * Varu[gen, t]
    )
    @constraint(
        model,
        ConstrRampUp[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t] - Varp[gen, t-1] <=
        RampUp[gen] * Varu[gen, t-1] + StartUp[gen] * Varv[gen, t]
    )
    @constraint(
        model,
        ConstrRampDown[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t-1] - Varp[gen, t] <=
        RampDown[gen] * Varu[gen, t] + ShutDown[gen] * Varw[gen, t]
    )

    @objective(
        model,
        Min,
        sum(
            sum(
                NoLoadConsumption[gen] * Varu[gen, t] +
                FixedCost[gen] * Varv[gen, t] +
                MarginalCost[gen] * Varp[gen, t] for t = 1:T
            ) for gen = 1:NbGen
        ) + LostLoad * sum(L[t] - VarL[t] for t = 1:T)
    )
    optimize!(model)

    prices = shadow_price.(loads) # prices = ProjBox(dual.(loads), PCD, PCU)
    return prices
end

function GetShift(instance)
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

    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model)

    @variable(model, Varp[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varpbar[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varu[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)
    @variable(model, Varv[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)
    @variable(model, Varw[g = 1:NbGen, t = 0:T+1], lower_bound = 0, upper_bound = 1)

    @variable(model, VarL[t = 1:T], lower_bound = 0)
    @constraint(model, [t = 1:T], VarL[t] <= L[t])

    @constraint(model, loads[t = 1:T], sum(Varp[g, t] for g = 1:NbGen) == VarL[t]) # Balance constraint

    @constraint(
        model,
        ConstrLogical[gen = 1:NbGen, t = 1:T+1],
        Varu[gen, t] - Varu[gen, t-1] == Varv[gen, t] - Varw[gen, t]
    )
    @constraint(
        model,
        ConstrMinUpTime[gen = 1:NbGen, t = UpTime[gen]:T+1],
        sum(Varv[gen, i] for i = t-UpTime[gen]+1:t) <= Varu[gen, t]
    )
    @constraint(
        model,
        ConstrMinDownTime[gen = 1:NbGen, t = DownTime[gen]:T+1],
        sum(Varw[gen, i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[gen, t])
    )
    @constraint(
        model,
        ConstrGenLimits1[gen = 1:NbGen, t = 0:T+1],
        MinRunCapacity[gen] * Varu[gen, t] <= Varp[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits2[gen = 1:NbGen, t = 0:T+1],
        Varp[gen, t] <= Varpbar[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits3[gen = 1:NbGen, t = 0:T+1],
        Varpbar[gen, t] <= MaxRunCapacity[gen] * Varu[gen, t]
    )
    @constraint(
        model,
        ConstrRampUp[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t] - Varp[gen, t-1] <=
        RampUp[gen] * Varu[gen, t-1] + StartUp[gen] * Varv[gen, t]
    )
    @constraint(
        model,
        ConstrRampDown[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t-1] - Varp[gen, t] <=
        RampDown[gen] * Varu[gen, t] + ShutDown[gen] * Varw[gen, t]
    )

    @objective(
        model,
        Min,
        sum(
            sum(
                NoLoadConsumption[gen] * Varu[gen, t] +
                FixedCost[gen] * Varv[gen, t] +
                MarginalCost[gen] * Varp[gen, t] for t = 1:T
            ) for gen = 1:NbGen
        ) + LostLoad * (sum(L[t] - VarL[t] for t = 1:T))
    )
    optimize!(model)

    U = value.(model[:Varu]).data
    V = value.(model[:Varv]).data
    W = value.(model[:Varw]).data
    P = value.(model[:Varp]).data
    Pbar = value.(model[:Varpbar]).data
    LL = value.(model[:VarL])
    return U, V, W, P, Pbar, LL
end

function LagrangianFunction(instance, prices)
    # TODO: Add ref to the name of each constraint in the paper

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

    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model)

    @variable(model, Varp[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varpbar[g = 1:NbGen, t = 0:T+1], lower_bound = 0)
    @variable(model, Varu[g = 1:NbGen, t = 0:T+1], Bin)
    @variable(model, Varv[g = 1:NbGen, t = 0:T+1], Bin)
    @variable(model, Varw[g = 1:NbGen, t = 0:T+1], Bin)

    @variable(model, VarL[t = 1:T], lower_bound = 0)
    @constraint(model, [t = 1:T], VarL[t] <= L[t])

    @constraint(
        model,
        ConstrLogical[gen = 1:NbGen, t = 1:T+1],
        Varu[gen, t] - Varu[gen, t-1] == Varv[gen, t] - Varw[gen, t]
    )
    @constraint(
        model,
        ConstrMinUpTime[gen = 1:NbGen, t = UpTime[gen]:T+1],
        sum(Varv[gen, i] for i = t-UpTime[gen]+1:t) <= Varu[gen, t]
    )
    @constraint(
        model,
        ConstrMinDownTime[gen = 1:NbGen, t = DownTime[gen]:T+1],
        sum(Varw[gen, i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[gen, t])
    )
    @constraint(
        model,
        ConstrGenLimits1[gen = 1:NbGen, t = 0:T+1],
        MinRunCapacity[gen] * Varu[gen, t] <= Varp[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits2[gen = 1:NbGen, t = 0:T+1],
        Varp[gen, t] <= Varpbar[gen, t]
    )
    @constraint(
        model,
        ConstrGenLimits3[gen = 1:NbGen, t = 0:T+1],
        Varpbar[gen, t] <= MaxRunCapacity[gen] * Varu[gen, t]
    )
    @constraint(
        model,
        ConstrRampUp[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t] - Varp[gen, t-1] <=
        RampUp[gen] * Varu[gen, t-1] + StartUp[gen] * Varv[gen, t]
    )
    @constraint(
        model,
        ConstrRampDown[gen = 1:NbGen, t = 1:T+1],
        Varpbar[gen, t-1] - Varp[gen, t] <=
        RampDown[gen] * Varu[gen, t] + ShutDown[gen] * Varw[gen, t]
    )

    @objective(
        model,
        Min,
        sum(
            sum(
                NoLoadConsumption[gen] * Varu[gen, t] +
                FixedCost[gen] * Varv[gen, t] +
                MarginalCost[gen] * Varp[gen, t] for t = 1:T
            ) for gen = 1:NbGen
        ) + LostLoad * (L[t] - sum(VarL[t] for t = 1:T)) -
        sum(prices[t] * (sum(Varp[gen, t] for gen = 1:NbGen) - VarL[t]) for t = 1:T)
    )
    optimize!(model)

    Varu, Varv, Varw, Varp, Varpbar, VarL = value.(Varu),
    value.(Varv),
    value.(Varw),
    value.(Varp),
    value.(Varpbar),
    value.(VarL)
    return LagrangianSolution(
        objective_value(model),
        Varu.data,
        Varv.data,
        Varw.data,
        Varp.data,
        Varpbar.data,
        VarL,
    )
end

function exact_oracle(instance, prices)
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

    ObjOracle = 0
    GradOracle = zeros(T)
    solutions = []
    my_lock = Threads.ReentrantLock()
    Threads.@threads for gen = 1:NbGen
        model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
        set_silent(model)
        set_optimizer_attributes(model, "MIPGap" => 1e-8, "MIPGapAbs" => 0)

        @variable(model, Varp[t = 0:T+1], lower_bound = 0)
        @variable(model, Varpbar[t = 0:T+1], lower_bound = 0)
        @variable(model, Varu[t = 0:T+1], Bin)
        @variable(model, Varv[t = 0:T+1], Bin)
        @variable(model, Varw[t = 0:T+1], Bin)

        @variable(model, VarL[t = 1:T], lower_bound = 0)
        @constraint(model, [t = 1:T], VarL[t] <= L[t])

        @constraint(
            model,
            ConstrLogical[t = 1:T+1],
            Varu[t] - Varu[t-1] == Varv[t] - Varw[t]
        )
        @constraint(model, ConstrGenLimits2[t = 0:T+1], Varp[t] <= Varpbar[t])

        @constraint(
            model,
            ConstrMinUpTime[t = UpTime[gen]:T+1],
            sum(Varv[i] for i = t-UpTime[gen]+1:t) <= Varu[t]
        )
        @constraint(
            model,
            ConstrMinDownTime[t = DownTime[gen]:T+1],
            sum(Varw[i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[t])
        )
        @constraint(
            model,
            ConstrGenLimits1[t = 0:T+1],
            MinRunCapacity[gen] * Varu[t] <= Varp[t]
        )
        @constraint(
            model,
            ConstrGenLimits3[t = 0:T+1],
            Varpbar[t] <= MaxRunCapacity[gen] * Varu[t]
        )
        @constraint(
            model,
            ConstrRampUp[t = 1:T+1],
            Varpbar[t] - Varp[t-1] <= RampUp[gen] * Varu[t-1] + StartUp[gen] * Varv[t]
        )
        @constraint(
            model,
            ConstrRampDown[t = 1:T+1],
            Varpbar[t-1] - Varp[t] <= RampDown[gen] * Varu[t] + ShutDown[gen] * Varw[t]
        )
        @objective(
            model,
            Min,
            sum(
                NoLoadConsumption[gen] * Varu[t] +
                FixedCost[gen] * Varv[t] +
                MarginalCost[gen] * Varp[t] +
                (LostLoad / NbGen) * (L[t] - VarL[t]) -
                prices[t] * (Varp[t] - (VarL[t] / NbGen)) for t = 1:T
            )
        )
        optimize!(model)

        ValL = value.(model[:VarL])
        Valp = value.(model[:Varp])
        Threads.lock(my_lock) do
            ObjOracle += objective_value(model)
            GradOracle += Array([ValL[t] / NbGen - Valp[t] for t = 1:T])
        end
    end
    return ObjOracle, GradOracle
end

function stochastic_oracle(instance, prices, batchsize)
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

    ObjOracle = 0
    GradOracle = zeros(T)
    indexes = randperm(NbGen)[1:batchsize]
    generators = [k for k = 1:NbGen]
    for gen in generators[indexes]
        model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
        set_silent(model)
        set_optimizer_attributes(model, "MIPGap" => 0, "MIPGapAbs" => 1e-8)

        @variable(model, Varp[t = 0:T+1], lower_bound = 0)
        @variable(model, Varpbar[t = 0:T+1], lower_bound = 0)
        @variable(model, Varu[t = 0:T+1], Bin)
        @variable(model, Varv[t = 0:T+1], Bin)
        @variable(model, Varw[t = 0:T+1], Bin)

        @variable(model, VarL[t = 1:T], lower_bound = 0)
        @constraint(model, [t = 1:T], VarL[t] <= L[t])

        @constraint(
            model,
            ConstrLogical[t = 1:T+1],
            Varu[t] - Varu[t-1] == Varv[t] - Varw[t]
        )
        @constraint(model, ConstrGenLimits2[t = 0:T+1], Varp[t] <= Varpbar[t])

        @constraint(
            model,
            ConstrMinUpTime[t = UpTime[gen]:T+1],
            sum(Varv[i] for i = t-UpTime[gen]+1:t) <= Varu[t]
        )
        @constraint(
            model,
            ConstrMinDownTime[t = DownTime[gen]:T+1],
            sum(Varw[i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[t])
        )
        @constraint(
            model,
            ConstrGenLimits1[t = 0:T+1],
            MinRunCapacity[gen] * Varu[t] <= Varp[t]
        )
        @constraint(
            model,
            ConstrGenLimits3[t = 0:T+1],
            Varpbar[t] <= MaxRunCapacity[gen] * Varu[t]
        )
        @constraint(
            model,
            ConstrRampUp[t = 1:T+1],
            Varpbar[t] - Varp[t-1] <= RampUp[gen] * Varu[t-1] + StartUp[gen] * Varv[t]
        )
        @constraint(
            model,
            ConstrRampDown[t = 1:T+1],
            Varpbar[t-1] - Varp[t] <= RampDown[gen] * Varu[t] + ShutDown[gen] * Varw[t]
        )

        @objective(
            model,
            Min,
            sum(
                NoLoadConsumption[gen] * Varu[t] +
                FixedCost[gen] * Varv[t] +
                MarginalCost[gen] * Varp[t] +
                (LostLoad / NbGen) * (L[t] - VarL[t]) -
                prices[t] * (Varp[t] - (VarL[t] / NbGen)) for t = 1:T
            )
        )
        optimize!(model)

        ObjOracle += objective_value(model)
        ValL = value.(model[:VarL])
        Valp = value.(model[:Varp])
        GradOracle += Array([ValL[t] / NbGen - Valp[t] for t = 1:T])
    end
    return ObjOracle, GradOracle
end

function exact_smooth_oracle(instance, prices, smoothing_parameter)
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

    ObjOracle = 0
    GradOracle = zeros(T)
    A = [Matrix(1.0I, T, T); zeros(Float64, 4 * T, T); Matrix((-1 / NbGen)I, T, T)]
    for gen = 1:NbGen
        model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
        set_silent(model)
        set_optimizer_attributes(model, "MIPGap" => 0, "MIPGapAbs" => 1e-8)

        @variable(model, Varp[t = 0:T+1], lower_bound = 0)
        @variable(model, Varpbar[t = 0:T+1], lower_bound = 0)
        @variable(model, Varu[t = 0:T+1], Bin)
        @variable(model, Varv[t = 0:T+1], Bin)
        @variable(model, Varw[t = 0:T+1], Bin)

        @variable(model, VarL[t = 1:T], lower_bound = 0)
        @constraint(model, [t = 1:T], VarL[t] <= L[t])

        @constraint(
            model,
            ConstrLogical[t = 1:T+1],
            Varu[t] - Varu[t-1] == Varv[t] - Varw[t]
        )
        @constraint(model, ConstrGenLimits2[t = 0:T+1], Varp[t] <= Varpbar[t])

        @constraint(
            model,
            ConstrMinUpTime[t = UpTime[gen]:T+1],
            sum(Varv[i] for i = t-UpTime[gen]+1:t) <= Varu[t]
        )
        @constraint(
            model,
            ConstrMinDownTime[t = DownTime[gen]:T+1],
            sum(Varw[i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[t])
        )
        @constraint(
            model,
            ConstrGenLimits1[t = 0:T+1],
            MinRunCapacity[gen] * Varu[t] <= Varp[t]
        )
        @constraint(
            model,
            ConstrGenLimits3[t = 0:T+1],
            Varpbar[t] <= MaxRunCapacity[gen] * Varu[t]
        )
        @constraint(
            model,
            ConstrRampUp[t = 1:T+1],
            Varpbar[t] - Varp[t-1] <= RampUp[gen] * Varu[t-1] + StartUp[gen] * Varv[t]
        )
        @constraint(
            model,
            ConstrRampDown[t = 1:T+1],
            Varpbar[t-1] - Varp[t] <= RampDown[gen] * Varu[t] + ShutDown[gen] * Varw[t]
        )

        prox = sum(
            Varp[t]^2 + Varpbar[t]^2 + Varu[t]^2 + Varv[t]^2 + Varw[t]^2 + VarL[t]^2 for
            t = 1:T
        )
        @objective(
            model,
            Min,
            sum(
                NoLoadConsumption[gen] * Varu[t] +
                FixedCost[gen] * Varv[t] +
                MarginalCost[gen] * Varp[t] +
                (LostLoad / NbGen) * (L[t] - VarL[t]) -
                prices[t] * (Varp[t] - (VarL[t] / NbGen)) +
                (smoothing_parameter / 2) * prox for t = 1:T
            )
        )
        optimize!(model)

        ObjOracle += objective_value(model)
        ValL = value.(model[:VarL])[1:T]
        Valp = value.(model[:Varp]).data[2:T+1]
        Valpbar = value.(model[:Varpbar]).data[2:T+1]
        Valu = value.(model[:Varu]).data[2:T+1]
        Valv = value.(model[:Varv]).data[2:T+1]
        Valw = value.(model[:Varw]).data[2:T+1]
        U = vcat(Valp, Valpbar, Valu, Valv, Valw, ValL)
        GradOracle += transpose(A) * U# Array([ValL[t]/NbGen - Valp[t] for t=1:T])
    end
    return ObjOracle, -GradOracle
end

function exact_translate_smooth_oracle(instance, prices, smoothing_parameter, shift)
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
    shiftU, shiftV, shiftW, shiftP, shiftPbar, shiftL = shift
    ObjOracle = 0
    GradOracle = zeros(T)
    A = [Matrix(1.0I, T, T); zeros(Float64, 4 * T, T); Matrix((-1 / NbGen)I, T, T)]
    for gen = 1:NbGen
        model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
        set_silent(model)
        set_optimizer_attributes(model, "MIPGap" => 0, "MIPGapAbs" => 1e-8)

        @variable(model, Varp[t = 0:T+1], lower_bound = 0)
        @variable(model, Varpbar[t = 0:T+1], lower_bound = 0)
        @variable(model, Varu[t = 0:T+1], Bin)
        @variable(model, Varv[t = 0:T+1], Bin)
        @variable(model, Varw[t = 0:T+1], Bin)

        @variable(model, VarL[t = 1:T], lower_bound = 0)
        @constraint(model, [t = 1:T], VarL[t] <= L[t])

        @constraint(
            model,
            ConstrLogical[t = 1:T+1],
            Varu[t] - Varu[t-1] == Varv[t] - Varw[t]
        )
        @constraint(model, ConstrGenLimits2[t = 0:T+1], Varp[t] <= Varpbar[t])

        @constraint(
            model,
            ConstrMinUpTime[t = UpTime[gen]:T+1],
            sum(Varv[i] for i = t-UpTime[gen]+1:t) <= Varu[t]
        )
        @constraint(
            model,
            ConstrMinDownTime[t = DownTime[gen]:T+1],
            sum(Varw[i] for i = t-DownTime[gen]+1:t) <= (1 - Varu[t])
        )
        @constraint(
            model,
            ConstrGenLimits1[t = 0:T+1],
            MinRunCapacity[gen] * Varu[t] <= Varp[t]
        )
        @constraint(
            model,
            ConstrGenLimits3[t = 0:T+1],
            Varpbar[t] <= MaxRunCapacity[gen] * Varu[t]
        )
        @constraint(
            model,
            ConstrRampUp[t = 1:T+1],
            Varpbar[t] - Varp[t-1] <= RampUp[gen] * Varu[t-1] + StartUp[gen] * Varv[t]
        )
        @constraint(
            model,
            ConstrRampDown[t = 1:T+1],
            Varpbar[t-1] - Varp[t] <= RampDown[gen] * Varu[t] + ShutDown[gen] * Varw[t]
        )

        prox = sum(
            (Varp[t] - shiftP[t])^2 +
            (Varpbar[t] - shiftPbar[t])^2 +
            (Varu[t] - shiftU[t])^2 +
            (Varv[t] - shiftV[t])^2 +
            (Varw[t] - shiftW[t])^2 +
            (VarL[t] - shiftL[t])^2 for t = 1:T
        )
        @objective(
            model,
            Min,
            sum(
                NoLoadConsumption[gen] * Varu[t] +
                FixedCost[gen] * Varv[t] +
                MarginalCost[gen] * Varp[t] +
                (LostLoad / NbGen) * (L[t] - VarL[t]) -
                prices[t] * (Varp[t] - (VarL[t] / NbGen)) +
                (smoothing_parameter / 2) * prox for t = 1:T
            )
        )
        optimize!(model)

        ObjOracle += objective_value(model)
        ValL = value.(model[:VarL])[1:T]
        Valp = value.(model[:Varp]).data[2:T+1]
        Valpbar = value.(model[:Varpbar]).data[2:T+1]
        Valu = value.(model[:Varu]).data[2:T+1]
        Valv = value.(model[:Varv]).data[2:T+1]
        Valw = value.(model[:Varw]).data[2:T+1]
        U = vcat(Valp, Valpbar, Valu, Valv, Valw, ValL)
        GradOracle += transpose(A) * U# Array([ValL[t]/NbGen - Valp[t] for t=1:T])
    end
    return ObjOracle, -GradOracle
end
