BLM = Dict("BE" => .99, "CA" => .98)
BPM = Dict("BE" => (.98, 1e4), "CA" => (.98, 5e4))
DA = Dict("BE" => 5e-2, "CA" => 5e-1) 
DOWG = Dict("BE" => 1e-4, "CA" => 1e-4) # CA is not optimized 
EPOL = Dict("BE" => 1.0, "CA" => 1e4) 
GD = Dict("BE" => 5e-4, "CA" => 5e-4) 
SGD = Dict("BE" => 1e-3, "CA" => 5e-4) # for batch of size 10
OGM = Dict("BE" => 1, "CA" => 2) # WTF ? Does not end.
CG = Dict("BE" => nothing, "CA" => nothing) # No hyperparameters
CRG = Dict("BE" => nothing, "CA" => nothing) # No hyperparameters
