Base.@kwdef mutable struct ThermalGen
    MinRunCapacity;
    MaxRunCapacity;
    RampUp;
    RampDown;
    StartUp;
    ShutDown;
    UpTime;
    DownTime;
    NoLoadConsumption;
    MarginalCost;
    FixedCost;
end

Base.@kwdef mutable struct Instance
    LostLoad;
    Load;
    ThermalGen::ThermalGen;
end

Base.@kwdef struct Solution
    Instance;
    ObjVal;
    ObjMatch;
    Prices;
    TotalRuntime;
    SolvingTime;
    Varp
    Varu
    Varv
    Varl
end

Base.@kwdef struct LagrangianSolution
    Obj;
    Varu;
    Varv;
    Varw;
    Varp;
    Varpbar;
    VarL;
end

Base.@kwdef struct MatchingSolution
    Obj;
    Varu;
    Varv;
    Varw;
    Varp;
    Varpbar;
    VarL;
end

struct SubProblem # used for ColumnGeneration
    model;
    Varp;
    Varpbar;
    Varu;
    Varv;
    Varw;
    VarCost;
    ConstrLogical;
    ConstrMinUpTime;
    ConstrMinDownTime;
    ConstrGenLimits1;
    ConstrGenLimits2;
    ConstrGenLimits3;
    ConstrRampUp;
    ConstrRampDown;
end

struct RestrictedMasterProgram
    model;
    VarZ;
    VarL;
    ConstrBalance;
    ConstrConvexComb;
end