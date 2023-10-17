# Metrics for the practicioner

function duality_gap(solution)
    obj_dual = LagrangianFunction(solution.Instance, solution.Prices)
    return solution.Obj_match - obj_dual
end

function total_cost(solution)
    NbGen = length(solution.Instance.ThermalGen.MinRunCapacity)
    TotalCosts = zeros(NbGen)
    for generator=1:NbGen
        NLC = solution.Instance.ThermalGen.NoLoadConsumption[generator]
        MC = solution.Instance.ThermalGen.MarginalCost[generator]
        FC = solution.Instance.ThermalGen.FixedCost[generator]
        Varu = solution.Varu[generator]
        Varv = solution.Varv[generator]
        Varp = solution.Varp[generator]
        T = length(solution.instance.Load)
        TotalCosts[generator] = sum((NLC * MC * Varu[t + 1] + FC * Varv[t + 1] + MC * Varp[t + 1]) for t=1:T)
    end
    return TotalCosts
end

function revenues(solution)
    NbGen = length(solution.Instance.ThermalGen.MinRunCapacity)
    Revenues = zeros(NbGen)
    for generator=1:NbGen
        T = length(solution.instance.Load)
        Revenues[generator] = sum(solution.Prices[t] * solution.Varp[generator, t + 1] for t=1:T)
    end
    return Revenues
end

function profits(solution)
    return revenues(solution) - total_cost(solution)
end

function selfish_profits(solution)
    NbGen = length(solution.Instance.ThermalGen.MinRunCapacity)
    SelfishProfits = zeros(NbGen)
    for generator=1:NbGen
        MaxProfit = MaxProfit_Producer(solution.Instance, solution.Prices, generator)
        SelfishProfits[generator] = MaxProfit
    end
    return SelfishProfits
end

function uplifts(solution)
    T = length(solution.instance.Load)
    MaxProfitConsumer = MaxProfit_Consumer(solution.Instance, solution.Prices)
    ReceivedProfitConsumer = sum((instance.LostLoad - solution.prices[t]) * solution.Varl[t] for t=1:T)
    UpliftConsumer = MaxProfitConsumer - ReceivedProfitConsumer
    UpliftProducer = selfish_profits(solution) - profits(solution)
    return vcat(UpliftProducer, UpliftConsumer)
end