using JuMP, SparseArrays

function set_T(UT, T_max, NbGen)
    index_gen = zeros(NbGen)
    for g = 1:NbGen
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
    A = sparse(zeros(NbGen, convert(Int64, maximum(index_gen))))
    B = sparse(zeros(NbGen, convert(Int64, maximum(index_gen))))

    for g = 1:NbGen
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

function Compute_A_B(U, NbGen)
    # Initialize an array to store the number of intervals per generation
    NbIntervalsGen = zeros(Int64, NbGen)

    # Check if there is more than one generator
    if NbGen > 1
        # Loop through each generator
        for g = 1:NbGen
            flag = false
            # Loop through each U value of the current generator
            for elem in U[g, :]
                if elem > 0
                    if flag == false
                        flag = true
                        NbIntervalsGen[g] += 1
                    end
                elseif flag == true
                    flag = false
                end
            end
        end
        # same logic without the outer loop
    else
        flag = false
        for elem in U
            if elem > 0
                if flag == false
                    flag = true
                    NbIntervalsGen[1] += 1
                end
            elseif flag == true
                flag = false
            end
        end
    end
    # Initialize arrays A and B to store interval start and end indices
    startIndexes = fill(-1, NbGen, maximum(NbIntervalsGen))
    stopIndexes = fill(-1, NbGen, maximum(NbIntervalsGen))
    a = 0
    if NbGen > 1
        for g = 1:NbGen
            flag = false
            count = 0
            index = 1
            for elem in U[g, :]
                count += 1
                if elem > 0
                    if flag == false
                        flag = true
                        a = count
                    end
                elseif flag == true
                    stopIndexes[g, index] = count - 1
                    startIndexes[g, index] = a
                    index += 1
                    flag = false
                end
            end
        end
        # Same logic without the outer loop
    else
        flag = false
        count = 0
        index = 1
        for elem in U
            count += 1
            if elem > 0
                if flag == false
                    flag = true
                    a = count
                end
            elseif flag == true
                stopIndexes[1, index] = count - 1
                startIndexes[1, index] = a
                index += 1
                flag = false
            end
        end
    end
    return startIndexes, stopIndexes, NbIntervalsGen
end


function A_B_In_Intervals(a, b, g, A, B, NbIntervalsGen)
    for i = 1:NbIntervalsGen[g]
        if a == A[g, i] && b == B[g, i]
            return true
        end
    end
    return false
end

function Update_A_B(U, NbGen, A, B, NbIntervalsGen)
    a = 0
    b = 0
    addedNbIntervals = zeros(Int64, NbGen)
    oldNbIntervalsGen = copy(NbIntervalsGen)
    if NbGen > 1
        for g = 1:NbGen
            flag = false
            count = 0
            for elem in U[g, :]
                count += 1
                if elem > if flag == false
                    flag = true
                    a = count
                end
                elseif flag == true
                    b = count - 1
                    if A_B_In_Intervals(a, b, g, A, B, oldNbIntervalsGen) == false
                        NbIntervalsGen[g] += 1
                        addedNbIntervals[g] += 1
                    end
                    flag = false
                end
            end
        end
    else
        flag = false
        count = 0
        for elem in U
            count += 1
            if elem > 0
                if flag == false
                    flag = true
                    a = count
                end
            elseif flag == true
                b = count - 1
                if A_B_In_Intervals(a, b, 1, A, B, oldNbIntervalsGen) == false
                    NbIntervalsGen[1] += 1
                    addedNbIntervals[1] += 1
                end
                flag = false
            end
        end
    end
    a = 0
    b = 0
    newA = (1) * ones(Int64, NbGen, maximum(NbIntervalsGen))
    newB = (1) * ones(Int64, NbGen, maximum(NbIntervalsGen))
    addedA = (1) * ones(Int64, NbGen, maximum(addedNbIntervals))
    addedB = (1) * ones(Int64, NbGen, maximum(addedNbIntervals))
    (i, j) = size(A)
    newA[1:i, 1:j] = A
    (i, j) = size(B)
    newB[1:i, 1:j] = B
    if NbGen > 1
        for g = 1:NbGen
            flag = false
            index = oldNbIntervalsGen[g] + 1
            index_added = 1
            count = 0
            for elem in U[g, :]
                count += 1
                if elem > 0
                    if flag == false
                        flag = true
                        a = count
                    end
                elseif flag == true
                    b = count - 1
                    if A_B_In_Intervals(a, b, g, A, B, oldNbIntervalsGen) == false
                        newA[g, index] = a
                        newB[g, index] = b
                        addedA[g, index_added] = a
                        addedB[g, index_added] = b
                        index += 1
                        index_added += 1
                    end
                    flag = false
                end
            end
        end
    else
        flag = false
        index = oldNbIntervalsGen[1] + 1
        index_added = 1
        count = 0
        for elem in U
            count += 1
            if elem > 0
                if flag == false
                    flag = true
                    a = count
                end
            elseif flag == true
                b = count - 1
                if A_B_In_Intervals(a, b, 1, A, B, oldNbIntervalsGen) == false
                    newA[1, index] = a
                    newB[1, index] = b
                    addedA[1, index_added] = a
                    addedB[1, index_added] = b
                    index += 1
                    index_added += 1
                end
                flag = false
            end
        end
    end
    return newA, newB, NbIntervalsGen, addedA, addedB, addedNbIntervals
end

function ColumnAndRowGeneration(instance, niter, eps)
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

    # Solving Matching problem, get Intervals
    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(model)
    set_optimizer_attribute(model, "Threads", 8)
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
        ConstrLogical[g = 1:NbGen, t = 1:T+1],
        Varu[g, t] - Varu[g, t-1] == Varv[g, t] - Varw[g, t]
    )
    @constraint(
        model,
        ConstrMinUpTime[g = 1:NbGen, t = UpTime[g]:T+1],
        sum(Varv[g, i] for i = t-UpTime[g]+1:t) <= Varu[g, t]
    )
    @constraint(
        model,
        ConstrMinDownTime[g = 1:NbGen, t = DownTime[g]:T+1],
        sum(Varw[g, i] for i = t-DownTime[g]+1:t) <= (1 - Varu[g, t])
    )
    @constraint(
        model,
        ConstrGenLimits1[g = 1:NbGen, t = 0:T+1],
        MinRunCapacity[g] * Varu[g, t] <= Varp[g, t]
    )
    @constraint(
        model,
        ConstrGenLimits2[g = 1:NbGen, t = 0:T+1],
        Varp[g, t] <= Varpbar[g, t]
    )
    @constraint(
        model,
        ConstrGenLimits3[g = 1:NbGen, t = 0:T+1],
        Varpbar[g, t] <= MaxRunCapacity[g] * Varu[g, t]
    )
    @constraint(
        model,
        ConstrRampUp[g = 1:NbGen, t = 1:T+1],
        Varpbar[g, t] - Varp[g, t-1] <= RampUp[g] * Varu[g, t-1] + StartUp[g] * Varv[g, t]
    )
    @constraint(
        model,
        ConstrRampDown[g = 1:NbGen, t = 1:T+1],
        Varpbar[g, t-1] - Varp[g, t] <= RampDown[g] * Varu[g, t] + ShutDown[g] * Varw[g, t]
    )

    for g = 1:NbGen
        JuMP.fix(model[:Varu][g, 0], 0; force = true)
        JuMP.fix(model[:Varv][g, 0], 0; force = true)
        JuMP.fix(model[:Varw][g, 0], 0; force = true)
        JuMP.fix(model[:Varp][g, 0], 0; force = true)
        JuMP.fix(model[:Varpbar][g, 0], 0; force = true)
    end

    @objective(
        model,
        Min,
        sum(
            sum(
                NoLoadConsumption[g] * model[:Varu][g, t] +
                FixedCost[g] * model[:Varv][g, t] +
                MarginalCost[g] * model[:Varp][g, t] for t = 1:T
            ) for g = 1:NbGen
        ) + LostLoad * sum(L[t] - model[:VarL][t] for t = 1:T)
    )
    optimize!(model)
    MatchingU = value.(model[:Varu]).data
    delete(model, loads)
    if NbGen > 1
        MatchingU = MatchingU[:, 2:T+1]
    else
        MatchingU = MatchingU[2:T+1]
    end
    (A, B, NbIntervalsGen) = set_T(UpTime, T, NbGen)
    addedA = copy(A)
    addedB = copy(B)
    addedNbIntervals = copy(NbIntervalsGen)

    # Restricted model
    RestrictedModel = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    set_silent(RestrictedModel)
    set_optimizer_attribute(RestrictedModel, "Threads", 8)
    @variable(RestrictedModel, 0 <= Varp[g = 1:NbGen, i = 1:addedNbIntervals[g], t = 0:T+1])
    @variable(
        RestrictedModel,
        0 <= Varpbar[g = 1:NbGen, i = 1:addedNbIntervals[g], t = 0:T+1]
    )
    @variable(RestrictedModel, 0 <= Vargamma[g = 1:NbGen, i = 1:addedNbIntervals[g]])
    @variable(RestrictedModel, 0 <= VarpTime[g = 1:NbGen, t = 0:T+1])
    @variable(RestrictedModel, 0 <= VarpbarTime[g = 1:NbGen, t = 0:T+1])
    @variable(RestrictedModel, 0 <= Varu[g = 1:NbGen, t = 0:T+1])
    @variable(RestrictedModel, 0 <= Varv[g = 1:NbGen, t = 0:T+1])
    @variable(RestrictedModel, 0 <= Varw[g = 1:NbGen, t = 0:T+1])
    @variable(RestrictedModel, 0 <= VarL[t = 1:T])
    @constraint(RestrictedModel, [t = 1:T], VarL[t] <= L[t])
    @constraint(
        RestrictedModel,
        [g = 1:NbGen, i = 1:addedNbIntervals[g]],
        Vargamma[g, i] <= 1
    )

    DictGamma = Dict()
    DictP = Dict()
    DictPbar = Dict()

    for g = 1:NbGen
        for i = 1:addedNbIntervals[g]
            DictGamma[g, i] = RestrictedModel[:Vargamma][g, i]
            for t = 0:T+1
                DictP[g, i, t] = RestrictedModel[:Varp][g, i, t]
                DictPbar[g, i, t] = RestrictedModel[:Varpbar][g, i, t]
            end
        end
    end

    CountP = 0
    CountGamma = 0
    CountPTime = 0
    CountU = 0
    for g = 1:NbGen
        CountPTime += T + 2
        CountU += T + 2
        for i = 1:addedNbIntervals[g]
            CountP += T + 2
            CountGamma += 1
        end
    end

    # Feasible Dispatch Polytope
    @constraint(
        RestrictedModel,
        no_production[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = 0:T+1;
            (t < addedA[g, i] || t > addedB[g, i]),
        ],
        DictP[g, i, t] <= 0
    )
    @constraint(
        RestrictedModel,
        min_output[g = 1:NbGen, i = 1:addedNbIntervals[g], t = addedA[g, i]:addedB[g, i]],
        -DictP[g, i, t] <= -DictGamma[g, i] * MinRunCapacity[g]
    )
    @constraint(
        RestrictedModel,
        max_output[g = 1:NbGen, i = 1:addedNbIntervals[g], t = addedA[g, i]:addedB[g, i]],
        DictP[g, i, t] - DictPbar[g, i, t] <= 0
    )
    @constraint(
        RestrictedModel,
        upper_bound_on_power_output_1[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]:addedB[g, i],
        ],
        DictPbar[g, i, t] <= DictGamma[g, i] * MaxRunCapacity[g]
    )
    @constraint(
        RestrictedModel,
        upper_bound_on_power_output_2[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]:addedB[g, i],
        ],
        DictPbar[g, i, t] <=
        DictGamma[g, i] * (StartUp[g] + (t - addedA[g, i]) * RampUp[g])
    )
    @constraint(
        RestrictedModel,
        upper_bound_on_power_output_3[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]:addedB[g, i],
        ],
        DictPbar[g, i, t] <=
        DictGamma[g, i] * (ShutDown[g] + (addedB[g, i] - t) * RampDown[g])
    )
    @constraint(
        RestrictedModel,
        limit_power_jumps_up_1[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]+1:addedB[g, i],
        ],
        DictPbar[g, i, t] - DictP[g, i, t-1] <= DictGamma[g, i] * RampUp[g]
    )
    @constraint(
        RestrictedModel,
        limit_power_jumps_up_2[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]+1:addedB[g, i],
        ],
        DictPbar[g, i, t] - DictP[g, i, t-1] <=
        DictGamma[g, i] *
        (ShutDown[g] + (addedB[g, i] - t) * RampDown[g] - MinRunCapacity[g])
    )
    @constraint(
        RestrictedModel,
        limit_power_jumps_down_1[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]+1:addedB[g, i],
        ],
        DictPbar[g, i, t-1] - DictP[g, i, t] <= DictGamma[g, i] * RampDown[g]
    )
    @constraint(
        RestrictedModel,
        limit_power_jumps_down_2[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]+1:addedB[g, i],
        ],
        DictPbar[g, i, t-1] - DictP[g, i, t] <=
        DictGamma[g, i] * (StartUp[g] + (t - addedA[g, i]) * RampUp[g] - MinRunCapacity[g])
    )
    @constraint(
        RestrictedModel,
        minimum_down_time_constraint[g = 1:NbGen, t = 1:T],
        sum(
            DictGamma[g, i] for i = 1:addedNbIntervals[g] if
            (addedA[g, i] <= t && t <= addedB[g, i] + DownTime[g])
        ) - 1 <= 0
    )
    @constraint(
        RestrictedModel,
        production_time[g = 1:NbGen, t = 1:T],
        VarpTime[g, t] - sum(DictP[g, i, t] for i = 1:addedNbIntervals[g]) == 0
    )
    @constraint(
        RestrictedModel,
        productionbar_time[g = 1:NbGen, t = 1:T],
        VarpbarTime[g, t] - sum(DictPbar[g, i, t] for i = 1:addedNbIntervals[g]) == 0
    )
    @constraint(
        RestrictedModel,
        commitment_status[g = 1:NbGen, t = 0:T+1],
        sum(
            DictGamma[g, i] for
            i = 1:addedNbIntervals[g] if (addedA[g, i] <= t && t <= addedB[g, i])
        ) - Varu[g, t] == 0
    )
    @constraint(
        RestrictedModel,
        startup_status[g = 1:NbGen, t = 0:T+1],
        sum(DictGamma[g, i] for i = 1:addedNbIntervals[g] if (t == addedA[g, i])) -
        Varv[g, t] == 0
    )
    @constraint(
        RestrictedModel,
        shutdown_status[g = 1:NbGen, t = 0:T+1],
        sum(DictGamma[g, i] for i = 1:addedNbIntervals[g] if (t == addedB[g, i] + 1)) -
        Varw[g, t] == 0
    )

    @constraint(
        RestrictedModel,
        loads[t = 1:T],
        sum(VarpTime[g, t] for g = 1:NbGen) - VarL[t] == 0
    )

    for g = 1:NbGen
        JuMP.fix(RestrictedModel[:Varu][g, 0], 0; force = true)
        JuMP.fix(RestrictedModel[:Varv][g, 0], 0; force = true)
        JuMP.fix(RestrictedModel[:Varw][g, 0], 0; force = true)
        JuMP.fix(RestrictedModel[:VarpTime][g, 0], 0; force = true)
        JuMP.fix(RestrictedModel[:VarpbarTime][g, 0], 0; force = true)
    end

    @objective(
        RestrictedModel,
        Min,
        sum(
            sum(
                NoLoadConsumption[g] * RestrictedModel[:Varu][g, t] +
                FixedCost[g] * RestrictedModel[:Varv][g, t] +
                MarginalCost[g] * RestrictedModel[:VarpTime][g, t] for t = 1:T
            ) for g = 1:NbGen
        ) + LostLoad * sum(L[t] - RestrictedModel[:VarL][t] for t = 1:T)
    )


    # Start iterations
    beta = 0
    arrObjRestricted = []
    arrObjPricing = []
    Prices = []

    for iter = 1:niter
        # Solve Restricted Problem
        optimize!(RestrictedModel)
        ObjRestricted = objective_value(RestrictedModel)
        Price = dual.(RestrictedModel[:loads])
        push!(Prices, Price)
        push!(arrObjRestricted, ObjRestricted)
        # Solve Pricing Problem
        @objective(
            model,
            Min,
            sum(
                sum(
                    NoLoadConsumption[g] * model[:Varu][g, t] +
                    FixedCost[g] * model[:Varv][g, t] +
                    MarginalCost[g] * model[:Varp][g, t] for t = 1:T
                ) for g = 1:NbGen
            ) + LostLoad * sum(L[t] - model[:VarL][t] for t = 1:T) - sum(
                Price[t] * (sum(model[:Varp][g, t] for g = 1:NbGen) - model[:VarL][t]) for
                t = 1:T
            )
        )
        optimize!(model)
        ObjPricing = objective_value(model)
        push!(arrObjPricing, ObjPricing)
        PricingU = value.(model[:Varu]).data

        if iter == 1
            beta = ObjPricing
        else
            beta = maximum([ObjPricing, beta])
        end

        if ObjRestricted <= beta + eps
            break
        end

        if NbGen > 1
            U = PricingU[:, 2:T+1]
        else
            U = PricingU[2:T+1]
        end

        newA, newB, newNbIntervals, addedA, addedB, addedNbIntervals =
            Update_A_B(U, NbGen, A, B, NbIntervalsGen)

        A = copy(newA)
        B = copy(newB)
        NbIntervalsGen = copy(newNbIntervals)

        # Update variables
        for g = 1:NbGen
            for i = 1:NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g]
                newVarp = @variable(RestrictedModel, [[g], [i], 0:T+1], lower_bound = 0)
                newVarpbar = @variable(RestrictedModel, [[g], [i], 0:T+1], lower_bound = 0)
                for t = 0:T+1
                    RestrictedModel[:Varp][g, i, t] = newVarp[g, i, t]
                    RestrictedModel[:Varpbar][g, i, t] = newVarpbar[g, i, t]
                end
                newVargamma = @variable(RestrictedModel, [[g], [i]], lower_bound = 0)
                RestrictedModel[:Vargamma][g, i] = newVargamma[g, i]
            end
        end

        # New constraints
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = 0:T+1;
                (t < A[g, i] || t > B[g, i]),
            ],
            RestrictedModel[:Varp][g, i, t] <= 0
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]:B[g, i],
            ],
            -RestrictedModel[:Varp][g, i, t] <=
            -RestrictedModel[:Vargamma][g, i] * MinRunCapacity[g]
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]:B[g, i],
            ],
            RestrictedModel[:Varp][g, i, t] - RestrictedModel[:Varpbar][g, i, t] <= 0
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]:B[g, i],
            ],
            RestrictedModel[:Varpbar][g, i, t] <=
            RestrictedModel[:Vargamma][g, i] * MaxRunCapacity[g]
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]:B[g, i],
            ],
            RestrictedModel[:Varpbar][g, i, t] <=
            RestrictedModel[:Vargamma][g, i] * (StartUp[g] + (t - A[g, i]) * RampUp[g])
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]:B[g, i],
            ],
            RestrictedModel[:Varpbar][g, i, t] <=
            RestrictedModel[:Vargamma][g, i] * (ShutDown[g] + (B[g, i] - t) * RampDown[g])
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]+1:B[g, i],
            ],
            RestrictedModel[:Varpbar][g, i, t] - RestrictedModel[:Varp][g, i, t-1] <=
            RestrictedModel[:Vargamma][g, i] * RampUp[g]
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]+1:B[g, i],
            ],
            RestrictedModel[:Varpbar][g, i, t] - RestrictedModel[:Varp][g, i, t-1] <=
            RestrictedModel[:Vargamma][g, i] *
            (ShutDown[g] + (B[g, i] - t) * RampDown[g] - MinRunCapacity[g])
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]+1:B[g, i],
            ],
            RestrictedModel[:Varpbar][g, i, t-1] - RestrictedModel[:Varp][g, i, t] <=
            RestrictedModel[:Vargamma][g, i] * RampDown[g]
        )
        @constraint(
            RestrictedModel,
            [
                g = 1:NbGen,
                i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                t = A[g, i]+1:B[g, i],
            ],
            RestrictedModel[:Varpbar][g, i, t-1] - RestrictedModel[:Varp][g, i, t] <=
            RestrictedModel[:Vargamma][g, i] *
            (StartUp[g] + (t - A[g, i]) * RampUp[g] - MinRunCapacity[g])
        )

        # Update constraints
        for g = 1:NbGen
            for t = 0:T+1
                for i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g]
                    if 1 <= t && t <= T
                        if A[g, i] <= t && t <= B[g, i] + DownTime[g]
                            set_normalized_coefficient(
                                RestrictedModel[:minimum_down_time_constraint][g, t],
                                RestrictedModel[:Vargamma][g, i],
                                1,
                            )
                        end
                        set_normalized_coefficient(
                            RestrictedModel[:production_time][g, t],
                            RestrictedModel[:Varp][g, i, t],
                            -1,
                        )
                        set_normalized_coefficient(
                            RestrictedModel[:productionbar_time][g, t],
                            RestrictedModel[:Varpbar][g, i, t],
                            -1,
                        )
                    end
                    if A[g, i] <= t && t <= B[g, i]
                        set_normalized_coefficient(
                            RestrictedModel[:commitment_status][g, t],
                            RestrictedModel[:Vargamma][g, i],
                            1,
                        )
                    end
                    if t == A[g, i]
                        set_normalized_coefficient(
                            RestrictedModel[:startup_status][g, t],
                            RestrictedModel[:Vargamma][g, i],
                            1,
                        )
                    end
                    if t == B[g, i] + 1
                        set_normalized_coefficient(
                            RestrictedModel[:shutdown_status][g, t],
                            RestrictedModel[:Vargamma][g, i],
                            1,
                        )
                    end
                end
            end
        end
    end

    return Prices[end], Prices, arrObjPricing
end

function tCRG(instance, budget, eps, verbose = -1)
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

    # Solving Matching problem, get Intervals
    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    # set_silent(model)
    set_optimizer_attribute(model, "Threads", 8)
    set_time_limit_sec(model, 25.0)
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
        ConstrLogical[g = 1:NbGen, t = 1:T+1],
        Varu[g, t] - Varu[g, t-1] == Varv[g, t] - Varw[g, t]
    )
    @constraint(
        model,
        ConstrMinUpTime[g = 1:NbGen, t = UpTime[g]:T+1],
        sum(Varv[g, i] for i = t-UpTime[g]+1:t) <= Varu[g, t]
    )
    @constraint(
        model,
        ConstrMinDownTime[g = 1:NbGen, t = DownTime[g]:T+1],
        sum(Varw[g, i] for i = t-DownTime[g]+1:t) <= (1 - Varu[g, t])
    )
    @constraint(
        model,
        ConstrGenLimits1[g = 1:NbGen, t = 0:T+1],
        MinRunCapacity[g] * Varu[g, t] <= Varp[g, t]
    )
    @constraint(
        model,
        ConstrGenLimits2[g = 1:NbGen, t = 0:T+1],
        Varp[g, t] <= Varpbar[g, t]
    )
    @constraint(
        model,
        ConstrGenLimits3[g = 1:NbGen, t = 0:T+1],
        Varpbar[g, t] <= MaxRunCapacity[g] * Varu[g, t]
    )
    @constraint(
        model,
        ConstrRampUp[g = 1:NbGen, t = 1:T+1],
        Varpbar[g, t] - Varp[g, t-1] <= RampUp[g] * Varu[g, t-1] + StartUp[g] * Varv[g, t]
    )
    @constraint(
        model,
        ConstrRampDown[g = 1:NbGen, t = 1:T+1],
        Varpbar[g, t-1] - Varp[g, t] <= RampDown[g] * Varu[g, t] + ShutDown[g] * Varw[g, t]
    )

    for g = 1:NbGen
        JuMP.fix(model[:Varu][g, 0], 0; force = true)
        JuMP.fix(model[:Varv][g, 0], 0; force = true)
        JuMP.fix(model[:Varw][g, 0], 0; force = true)
        JuMP.fix(model[:Varp][g, 0], 0; force = true)
        JuMP.fix(model[:Varpbar][g, 0], 0; force = true)
    end
    @objective(
        model,
        Min,
        sum(
            sum(
                NoLoadConsumption[g] * model[:Varu][g, t] +
                FixedCost[g] * model[:Varv][g, t] +
                MarginalCost[g] * model[:Varp][g, t] for t = 1:T
            ) for g = 1:NbGen
        ) + LostLoad * sum(L[t] - model[:VarL][t] for t = 1:T)
    )
    optimize!(model)
    MatchingU = value.(model[:Varu]).data
    delete(model, loads)
    if NbGen > 1
        MatchingU = MatchingU[:, 2:T+1]
    else
        MatchingU = MatchingU[2:T+1]
    end
    if verbose > 0
        @info " set T"
    end
    (A, B, NbIntervalsGen) = set_T(UpTime, T, NbGen) # Compute_A_B(MatchingU, NbGen) # 
    addedA = copy(A)
    addedB = copy(B)
    addedNbIntervals = copy(NbIntervalsGen)

    # Restricted model
    RestrictedModel = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV[]))
    # set_silent(RestrictedModel)
    set_optimizer_attribute(RestrictedModel, "Threads", 8)
    if verbose > 0
        @info "create restricted model"
    end
    @variable(RestrictedModel, 0 <= Varp[g = 1:NbGen, i = 1:addedNbIntervals[g], t = 0:T+1])
    @info "varp done"
    @variable(
        RestrictedModel,
        0 <= Varpbar[g = 1:NbGen, i = 1:addedNbIntervals[g], t = 0:T+1]
    )
    @info "varpbar done"
    @variable(RestrictedModel, 0 <= Vargamma[g = 1:NbGen, i = 1:addedNbIntervals[g]])
    @info "vargamma done"
    @variable(RestrictedModel, 0 <= VarpTime[g = 1:NbGen, t = 0:T+1])
    @info "varptime done"
    @variable(RestrictedModel, 0 <= VarpbarTime[g = 1:NbGen, t = 0:T+1])
    @info "varpbartime done"
    @variable(RestrictedModel, 0 <= Varu[g = 1:NbGen, t = 0:T+1])
    @info "varu done"
    @variable(RestrictedModel, 0 <= Varv[g = 1:NbGen, t = 0:T+1])
    @info "varv done"
    @variable(RestrictedModel, 0 <= Varw[g = 1:NbGen, t = 0:T+1])
    @info "varw done"
    @variable(RestrictedModel, 0 <= VarL[t = 1:T])
    @info "varL done"
    @constraint(RestrictedModel, [t = 1:T], VarL[t] <= L[t])
    @info "varL constraint done"
    @constraint(
        RestrictedModel,
        [g = 1:NbGen, i = 1:addedNbIntervals[g]],
        Vargamma[g, i] <= 1
    )
    @info "vargamma constraint done"
    DictGamma = Dict()
    DictP = Dict()
    DictPbar = Dict()
    if verbose > 0
        @info "Filling dicts."
    end
    for g = 1:NbGen
        for i = 1:addedNbIntervals[g]
            DictGamma[g, i] = Vargamma[g, i]
            for t = 0:T+1
                DictP[g, i, t] = Varp[g, i, t]
                DictPbar[g, i, t] = Varpbar[g, i, t]
            end
        end
    end
    CountP = 0
    CountGamma = 0
    CountPTime = 0
    CountU = 0
    if verbose > 0
        @info "filling counts"
    end
    for g = 1:NbGen
        CountPTime += T + 2
        CountU += T + 2
        for i = 1:addedNbIntervals[g]
            CountP += T + 2
            CountGamma += 1
        end
    end
    if verbose > 0
        @info "Feasible Dispatch Polytope"
    end
    # Feasible Dispatch Polytope
    @constraint(
        RestrictedModel,
        no_production[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = 0:T+1;
            (t < addedA[g, i] || t > addedB[g, i]),
        ],
        DictP[g, i, t] <= 0
    )
    @constraint(
        RestrictedModel,
        min_output[g = 1:NbGen, i = 1:addedNbIntervals[g], t = addedA[g, i]:addedB[g, i]],
        -DictP[g, i, t] <= -DictGamma[g, i] * MinRunCapacity[g]
    )
    @constraint(
        RestrictedModel,
        max_output[g = 1:NbGen, i = 1:addedNbIntervals[g], t = addedA[g, i]:addedB[g, i]],
        DictP[g, i, t] - DictPbar[g, i, t] <= 0
    )
    @constraint(
        RestrictedModel,
        upper_bound_on_power_output_1[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]:addedB[g, i],
        ],
        DictPbar[g, i, t] <= DictGamma[g, i] * MaxRunCapacity[g]
    )
    @constraint(
        RestrictedModel,
        upper_bound_on_power_output_2[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]:addedB[g, i],
        ],
        DictPbar[g, i, t] <=
        DictGamma[g, i] * (StartUp[g] + (t - addedA[g, i]) * RampUp[g])
    )
    @constraint(
        RestrictedModel,
        upper_bound_on_power_output_3[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]:addedB[g, i],
        ],
        DictPbar[g, i, t] <=
        DictGamma[g, i] * (ShutDown[g] + (addedB[g, i] - t) * RampDown[g])
    )
    @constraint(
        RestrictedModel,
        limit_power_jumps_up_1[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]+1:addedB[g, i],
        ],
        DictPbar[g, i, t] - DictP[g, i, t-1] <= DictGamma[g, i] * RampUp[g]
    )
    @constraint(
        RestrictedModel,
        limit_power_jumps_up_2[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]+1:addedB[g, i],
        ],
        DictPbar[g, i, t] - DictP[g, i, t-1] <=
        DictGamma[g, i] *
        (ShutDown[g] + (addedB[g, i] - t) * RampDown[g] - MinRunCapacity[g])
    )
    @constraint(
        RestrictedModel,
        limit_power_jumps_down_1[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]+1:addedB[g, i],
        ],
        DictPbar[g, i, t-1] - DictP[g, i, t] <= DictGamma[g, i] * RampDown[g]
    )
    @constraint(
        RestrictedModel,
        limit_power_jumps_down_2[
            g = 1:NbGen,
            i = 1:addedNbIntervals[g],
            t = addedA[g, i]+1:addedB[g, i],
        ],
        DictPbar[g, i, t-1] - DictP[g, i, t] <=
        DictGamma[g, i] * (StartUp[g] + (t - addedA[g, i]) * RampUp[g] - MinRunCapacity[g])
    )
    @constraint(
        RestrictedModel,
        minimum_down_time_constraint[g = 1:NbGen, t = 1:T],
        sum(
            DictGamma[g, i] for i = 1:addedNbIntervals[g] if
            (addedA[g, i] <= t && t <= addedB[g, i] + DownTime[g])
        ) - 1 <= 0
    )
    @constraint(
        RestrictedModel,
        production_time[g = 1:NbGen, t = 1:T],
        VarpTime[g, t] - sum(DictP[g, i, t] for i = 1:addedNbIntervals[g]) == 0
    )
    @constraint(
        RestrictedModel,
        productionbar_time[g = 1:NbGen, t = 1:T],
        VarpbarTime[g, t] - sum(DictPbar[g, i, t] for i = 1:addedNbIntervals[g]) == 0
    )
    @constraint(
        RestrictedModel,
        commitment_status[g = 1:NbGen, t = 0:T+1],
        sum(
            DictGamma[g, i] for
            i = 1:addedNbIntervals[g] if (addedA[g, i] <= t && t <= addedB[g, i])
        ) - Varu[g, t] == 0
    )
    @constraint(
        RestrictedModel,
        startup_status[g = 1:NbGen, t = 0:T+1],
        sum(DictGamma[g, i] for i = 1:addedNbIntervals[g] if (t == addedA[g, i])) -
        Varv[g, t] == 0
    )
    @constraint(
        RestrictedModel,
        shutdown_status[g = 1:NbGen, t = 0:T+1],
        sum(DictGamma[g, i] for i = 1:addedNbIntervals[g] if (t == addedB[g, i] + 1)) -
        Varw[g, t] == 0
    )

    @constraint(
        RestrictedModel,
        loads[t = 1:T],
        sum(VarpTime[g, t] for g = 1:NbGen) - VarL[t] == 0
    )
    if verbose > 0
        @info "fixing constraints"
    end
    for g = 1:NbGen
        JuMP.fix(RestrictedModel[:Varu][g, 0], 0; force = true)
        JuMP.fix(RestrictedModel[:Varv][g, 0], 0; force = true)
        JuMP.fix(RestrictedModel[:Varw][g, 0], 0; force = true)
        JuMP.fix(RestrictedModel[:VarpTime][g, 0], 0; force = true)
        JuMP.fix(RestrictedModel[:VarpbarTime][g, 0], 0; force = true)
    end
    if verbose > 0
        @info "dexfining objective"
    end
    @objective(
        RestrictedModel,
        Min,
        sum(
            sum(
                NoLoadConsumption[g] * RestrictedModel[:Varu][g, t] +
                FixedCost[g] * RestrictedModel[:Varv][g, t] +
                MarginalCost[g] * RestrictedModel[:VarpTime][g, t] for t = 1:T
            ) for g = 1:NbGen
        ) + LostLoad * sum(L[t] - RestrictedModel[:VarL][t] for t = 1:T)
    )


    # Start iterations
    beta = 0
    arrObjRestricted = []
    arrObjPricing = []
    Prices = []

    time_vector = [0.0]
    iter = 1
    if verbose > 0
        @info "Starting iterations"
    end
    while time_vector[end] <= budget
        @info "Iteration $iter"
        it_time = @elapsed begin
            # Solve Restricted Problem
            optimize!(RestrictedModel)
            ObjRestricted = objective_value(RestrictedModel)
            Price = dual.(RestrictedModel[:loads])
            push!(Prices, Price)
            push!(arrObjRestricted, ObjRestricted)
            # Solve Pricing Problem
            @objective(
                model,
                Min,
                sum(
                    sum(
                        NoLoadConsumption[g] * model[:Varu][g, t] +
                        FixedCost[g] * model[:Varv][g, t] +
                        MarginalCost[g] * model[:Varp][g, t] for t = 1:T
                    ) for g = 1:NbGen
                ) + LostLoad * sum(L[t] - model[:VarL][t] for t = 1:T) - sum(
                    Price[t] * (sum(model[:Varp][g, t] for g = 1:NbGen) - model[:VarL][t]) for
                    t = 1:T
                )
            )
            my_time = @elapsed begin
                optimize!(model)
            end
            if verbose > 0
                @info "Optimizing in iteration $iter took $my_time"
            end
            ObjPricing = objective_value(model)
            push!(arrObjPricing, ObjPricing)
            PricingU = value.(model[:Varu]).data

            if iter == 1
                beta = ObjPricing
            else
                beta = maximum([ObjPricing, beta])
            end

            if ObjRestricted <= beta + eps
                break
            end

            if NbGen > 1
                U = PricingU[:, 2:T+1]
            else
                U = PricingU[2:T+1]
            end

            newA, newB, newNbIntervals, addedA, addedB, addedNbIntervals =
                Update_A_B(U, NbGen, A, B, NbIntervalsGen)

            A = copy(newA)
            B = copy(newB)
            NbIntervalsGen = copy(newNbIntervals)

            # Update variables
            for g = 1:NbGen
                for i = 1:NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g]
                    newVarp = @variable(RestrictedModel, [[g], [i], 0:T+1], lower_bound = 0)
                    newVarpbar =
                        @variable(RestrictedModel, [[g], [i], 0:T+1], lower_bound = 0)
                    for t = 0:T+1
                        RestrictedModel[:Varp][g, i, t] = newVarp[g, i, t]
                        RestrictedModel[:Varpbar][g, i, t] = newVarpbar[g, i, t]
                    end
                    newVargamma =
                        @variable(RestrictedModel, [[g], [i]], lower_bound = 0)
                    RestrictedModel[:Vargamma][g, i] = newVargamma[g, i]
                end
            end

            # New constraints
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = 0:T+1;
                    (t < A[g, i] || t > B[g, i]),
                ],
                RestrictedModel[:Varp][g, i, t] <= 0
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]:B[g, i],
                ],
                -RestrictedModel[:Varp][g, i, t] <=
                -RestrictedModel[:Vargamma][g, i] * MinRunCapacity[g]
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]:B[g, i],
                ],
                RestrictedModel[:Varp][g, i, t] - RestrictedModel[:Varpbar][g, i, t] <= 0
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]:B[g, i],
                ],
                RestrictedModel[:Varpbar][g, i, t] <=
                RestrictedModel[:Vargamma][g, i] * MaxRunCapacity[g]
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]:B[g, i],
                ],
                RestrictedModel[:Varpbar][g, i, t] <=
                RestrictedModel[:Vargamma][g, i] * (StartUp[g] + (t - A[g, i]) * RampUp[g])
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]:B[g, i],
                ],
                RestrictedModel[:Varpbar][g, i, t] <=
                RestrictedModel[:Vargamma][g, i] *
                (ShutDown[g] + (B[g, i] - t) * RampDown[g])
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]+1:B[g, i],
                ],
                RestrictedModel[:Varpbar][g, i, t] - RestrictedModel[:Varp][g, i, t-1] <=
                RestrictedModel[:Vargamma][g, i] * RampUp[g]
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]+1:B[g, i],
                ],
                RestrictedModel[:Varpbar][g, i, t] - RestrictedModel[:Varp][g, i, t-1] <=
                RestrictedModel[:Vargamma][g, i] *
                (ShutDown[g] + (B[g, i] - t) * RampDown[g] - MinRunCapacity[g])
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]+1:B[g, i],
                ],
                RestrictedModel[:Varpbar][g, i, t-1] - RestrictedModel[:Varp][g, i, t] <=
                RestrictedModel[:Vargamma][g, i] * RampDown[g]
            )
            @constraint(
                RestrictedModel,
                [
                    g = 1:NbGen,
                    i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g],
                    t = A[g, i]+1:B[g, i],
                ],
                RestrictedModel[:Varpbar][g, i, t-1] - RestrictedModel[:Varp][g, i, t] <=
                RestrictedModel[:Vargamma][g, i] *
                (StartUp[g] + (t - A[g, i]) * RampUp[g] - MinRunCapacity[g])
            )

            # Update constraints
            for g = 1:NbGen
                for t = 0:T+1
                    for i = NbIntervalsGen[g]-addedNbIntervals[g]+1:NbIntervalsGen[g]
                        if 1 <= t && t <= T
                            if A[g, i] <= t && t <= B[g, i] + DownTime[g]
                                set_normalized_coefficient(
                                    RestrictedModel[:minimum_down_time_constraint][g, t],
                                    RestrictedModel[:Vargamma][g, i],
                                    1,
                                )
                            end
                            set_normalized_coefficient(
                                RestrictedModel[:production_time][g, t],
                                RestrictedModel[:Varp][g, i, t],
                                -1,
                            )
                            set_normalized_coefficient(
                                RestrictedModel[:productionbar_time][g, t],
                                RestrictedModel[:Varpbar][g, i, t],
                                -1,
                            )
                        end
                        if A[g, i] <= t && t <= B[g, i]
                            set_normalized_coefficient(
                                RestrictedModel[:commitment_status][g, t],
                                RestrictedModel[:Vargamma][g, i],
                                1,
                            )
                        end
                        if t == A[g, i]
                            set_normalized_coefficient(
                                RestrictedModel[:startup_status][g, t],
                                RestrictedModel[:Vargamma][g, i],
                                1,
                            )
                        end
                        if t == B[g, i] + 1
                            set_normalized_coefficient(
                                RestrictedModel[:shutdown_status][g, t],
                                RestrictedModel[:Vargamma][g, i],
                                1,
                            )
                        end
                    end
                end
            end
        end
        push!(time_vector, it_time + time_vector[end])
        iter += 1
    end
    f_iterates = Float64[]
    for ρ in Prices
        push!(f_iterates, Utilitaries.exact_oracle(instance, ρ)[1])
    end
    return Prices[end], Prices, f_iterates, time_vector
end
