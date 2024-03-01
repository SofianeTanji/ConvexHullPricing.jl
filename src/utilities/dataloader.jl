# Load JSON data as Instance object

function load_data(JSON_file)
    df = DataFrame(JSON3.read(JSON_file))
    Gen = df.thermal_generators[1]
    MinRunCapacity,
    MaxRunCapacity,
    RampUp,
    RampDown,
    UpTime,
    DownTime,
    StartUp,
    ShutDown,
    NoLoadConsumption,
    FixedCost,
    MarginalCost = [], [], [], [], [], [], [], [], [], [], []

    for (sym, gen) in Gen
        push!(MinRunCapacity, gen.power_output_minimum)
        push!(MaxRunCapacity, gen.power_output_maximum)
        push!(RampUp, gen.ramp_up_limit)
        push!(RampDown, gen.ramp_down_limit)
        push!(UpTime, gen.time_up_minimum)
        push!(DownTime, gen.time_down_minimum)

        if haskey(gen, "ramp_startup_limit")
            push!(StartUp, gen.ramp_startup_limit)
        else
            push!(StartUp, gen.power_output_minimum)
        end
        if haskey(gen, "ramp_shutdown_limit")
            push!(ShutDown, gen.ramp_shutdown_limit)
        else
            push!(ShutDown, gen.power_output_minimum)
        end
        if haskey(gen, "no_load_consumption")
            push!(NoLoadConsumption, gen.no_load_consumption)
        else
            push!(NoLoadConsumption, gen.time_down_minimum)
        end

        push!(FixedCost, gen.startup[1]["cost"])
        push!(MarginalCost, gen.piecewise_production[1]["cost"])
    end

    Gen = ThermalGen(
        MinRunCapacity = MinRunCapacity,
        MaxRunCapacity = MaxRunCapacity,
        RampUp = RampUp,
        RampDown = RampDown,
        StartUp = StartUp,
        ShutDown = ShutDown,
        UpTime = UpTime,
        DownTime = DownTime,
        NoLoadConsumption = NoLoadConsumption,
        MarginalCost = MarginalCost,
        FixedCost = FixedCost,
    )

    instance =
        Instance(LostLoad = 3000.0, Load = Array{Float64}(df.demand), ThermalGen = Gen)
    return instance
end

function load_rts_data(JSON_file)
    df = DataFrame(JSON3.read(JSON_file))
    Gen = df.thermal_generators[1]
    MinRunCapacity,
    MaxRunCapacity,
    RampUp,
    RampDown,
    UpTime,
    DownTime,
    StartUp,
    ShutDown,
    NoLoadConsumption,
    FixedCost,
    MarginalCost = [], [], [], [], [], [], [], [], [], [], []

    for (sym, gen) in Gen
        push!(MinRunCapacity, gen.power_output_minimum)
        push!(MaxRunCapacity, gen.power_output_maximum)
        push!(RampUp, gen.ramp_up_limit / 3.0)
        push!(RampDown, gen.ramp_down_limit / 3.0)
        push!(UpTime, gen.time_up_minimum)
        push!(DownTime, gen.time_down_minimum)

        if haskey(gen, "ramp_startup_limit")
            push!(StartUp, gen.ramp_startup_limit)
        else
            push!(StartUp, gen.power_output_minimum)
        end
        if haskey(gen, "ramp_shutdown_limit")
            push!(ShutDown, gen.ramp_shutdown_limit)
        else
            push!(ShutDown, gen.power_output_minimum)
        end
        if haskey(gen, "no_load_consumption")
            push!(NoLoadConsumption, gen.no_load_consumption)
        else
            push!(NoLoadConsumption, gen.time_down_minimum)
        end

        push!(FixedCost, gen.startup[1]["cost"])
        push!(MarginalCost, gen.piecewise_production[1]["cost"])
    end

    Gen = ThermalGen(
        MinRunCapacity = MinRunCapacity,
        MaxRunCapacity = MaxRunCapacity,
        RampUp = RampUp,
        RampDown = RampDown,
        StartUp = StartUp,
        ShutDown = ShutDown,
        UpTime = UpTime,
        DownTime = DownTime,
        NoLoadConsumption = NoLoadConsumption,
        MarginalCost = MarginalCost,
        FixedCost = FixedCost,
    )

    instance =
        Instance(LostLoad = 3000.0, Load = Array{Float64}(df.demand) ./ 2, ThermalGen = Gen)
    return instance
end

function load_ferc_data(JSON_file)
    df = DataFrame(JSON3.read(JSON_file))
    Gen = df.thermal_generators[1]
    MinRunCapacity,
    MaxRunCapacity,
    RampUp,
    RampDown,
    UpTime,
    DownTime,
    StartUp,
    ShutDown,
    NoLoadConsumption,
    FixedCost,
    MarginalCost = [], [], [], [], [], [], [], [], [], [], []

    for (sym, gen) in Gen
        push!(MinRunCapacity, gen.power_output_minimum)
        push!(MaxRunCapacity, gen.power_output_maximum)
        push!(RampUp, gen.ramp_up_limit)
        push!(RampDown, gen.ramp_down_limit)
        push!(UpTime, gen.time_up_minimum)
        push!(DownTime, gen.time_down_minimum)

        if haskey(gen, "ramp_startup_limit")
            push!(StartUp, gen.ramp_startup_limit)
        else
            push!(StartUp, gen.power_output_minimum)
        end
        if haskey(gen, "ramp_shutdown_limit")
            push!(ShutDown, gen.ramp_shutdown_limit)
        else
            push!(ShutDown, gen.power_output_minimum)
        end
        if haskey(gen, "no_load_consumption")
            push!(NoLoadConsumption, gen.no_load_consumption)
        else
            push!(NoLoadConsumption, gen.time_down_minimum)
        end

        push!(FixedCost, gen.startup[1]["cost"])
        push!(MarginalCost, gen.piecewise_production[1]["cost"])
    end

    Gen = ThermalGen(
        MinRunCapacity = MinRunCapacity,
        MaxRunCapacity = MaxRunCapacity,
        RampUp = RampUp,
        RampDown = RampDown,
        StartUp = StartUp,
        ShutDown = ShutDown,
        UpTime = UpTime,
        DownTime = DownTime,
        NoLoadConsumption = NoLoadConsumption,
        MarginalCost = MarginalCost,
        FixedCost = FixedCost,
    )

    instance = Instance(
        LostLoad = 3000.0,
        Load = Array{Float64}(df.demand) ./ 12,
        ThermalGen = Gen,
    )
    return instance
end
