cd(dirname(Base.source_path()))
using Distributions, Roots, QuadGK
include("./TVGG1PS_types.jl")
include("./TVGG1PS_functions_full_iter.jl")

# main
queue = "TVGG1PS"
control = "PD"
N = 10000

target_set = (0.1, 1.0, 10.0)
dist_set = ("ER","LN")
γ_set = (0.1, 0.01, 0.001)
T_set = (2000.0, 2000.0, 20000.0)

for dist in dist_set
    for target in target_set
        for i in 1:3
            arrival = dist
            service = dist
            coeff = (1.0, 0.2, γ_set[i])
            T = T_set[i]
            record = Record()
            do_experiment(queue, control, target, arrival, service, coeff, T, N, record)
        end
    end
end
