using Revise
using ConvexHullPricing
using DataFrames
using Plots
using JLD2
using JuMP, Gurobi
const UT = ConvexHullPricing.Utilitaries
const OPT = ConvexHullPricing.Optimizer

BEinstances = []
for file in readdir("data//belgian"; join=true)
  push!(BEinstances, UT.load_data(file))
end
instance = BEinstances[8]

XCRG, ITCRG, CRG, TCRG = OPT.tCRG(instance, 1., 1e-5, 1)