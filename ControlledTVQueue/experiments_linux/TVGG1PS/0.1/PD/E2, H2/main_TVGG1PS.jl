cd(dirname(Base.source_path()))
using Distributions, Roots, QuadGK
include("../../../types.jl")
include("../../../functions.jl")

# main
queue = "TVGG1PS"
target = 0.1
control = "PD"
arrival = "E2"
service = "H2"
N = 10000

γ_set = (0.1, 0.01, 0.001)
T_set = (200.0, 2000.0, 20000.0)
for i in 1:3
    coeff = (1.0, 0.2, γ_set[i])
    T = T_set[i]
    record = Record()
    do_experiment(queue, control, target, arrival, service, coeff, T, N, record)
end
