cd(dirname(Base.source_path()))
using Distributions

# main
queue = "GG1PS"
include("./$(queue)_types.jl")
include("./$(queue)_functions.jl")
arrival = "ER1"
service = "ER1"
service_rate = 1.0
T = 10000.0
N = 10000
record = Record()

do_experiment(queue, arrival, service, service_rate, T, N, record)
