Base.@kwdef mutable struct ThermalGen
    MinRunCapacity::Any
    MaxRunCapacity::Any
    RampUp::Any
    RampDown::Any
    StartUp::Any
    ShutDown::Any
    UpTime::Any
    DownTime::Any
    NoLoadConsumption::Any
    MarginalCost::Any
    FixedCost::Any
end

Base.@kwdef mutable struct Instance
    LostLoad::Any
    Load::Any
    ThermalGen::ThermalGen
end

Base.@kwdef struct Solution
    Instance::Any
    ObjVal::Any
    ObjMatch::Any
    Prices::Any
    TotalRuntime::Any
    SolvingTime::Any
    Varp::Any
    Varu::Any
    Varv::Any
    Varl::Any
end

Base.@kwdef struct LagrangianSolution
    Obj::Any
    Varu::Any
    Varv::Any
    Varw::Any
    Varp::Any
    Varpbar::Any
    VarL::Any
end

Base.@kwdef struct MatchingSolution
    Obj::Any
    Varu::Any
    Varv::Any
    Varw::Any
    Varp::Any
    Varpbar::Any
    VarL::Any
end

struct SubProblem # used for ColumnGeneration
    model::Any
    Varp::Any
    Varpbar::Any
    Varu::Any
    Varv::Any
    Varw::Any
    VarCost::Any
    ConstrLogical::Any
    ConstrMinUpTime::Any
    ConstrMinDownTime::Any
    ConstrGenLimits1::Any
    ConstrGenLimits2::Any
    ConstrGenLimits3::Any
    ConstrRampUp::Any
    ConstrRampDown::Any
end

struct RestrictedMasterProgram
    model::Any
    VarZ::Any
    VarL::Any
    ConstrBalance::Any
    ConstrConvexComb::Any
end
