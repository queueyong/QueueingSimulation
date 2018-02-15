cd(dirname(Base.source_path()))
using Distributions, Roots, QuadGK
include("./Types.jl")
include("./functions.jl")

# main
queue = "TVGG1PS"
control = "DM"
param = 0.25
arrival = "Exponential"
service = "Exponential"
coeff = (1.0, 0.2, 0.01)
T = 2000.0
N = 30
record = Record()
do_experiment(queue, control, param, arrival, service, coeff, T, N, record)
