cd(dirname(Base.source_path()))
using Distributions, Roots, QuadGK
include("./Types.jl")
include("./functions.jl")

# main
queue = "TVGG1"
control = "PSA"
param = 1.0
arrival = "Hyperexponential"
service = "Erlang"
coeff = (1.0, 0.2, 0.001)
T = 20000.0
N = 10
record = Record()
do_experiment(queue, control, param, arrival, service, coeff, T, N, record)
