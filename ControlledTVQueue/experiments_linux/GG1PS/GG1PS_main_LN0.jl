cd(dirname(Base.source_path()))
using Distributions

# main
queue = "GG1PS"
include("./$(queue)_types.jl")
include("./$(queue)_functions.jl")
arrival = "LN0"
service_rate_set = (1.0 , 2.0)

for service_rate in service_rate_set
    if service_rate == 1.0
        service_set = ("ER3","LN3")
    elseif service_rate == 2.0
        service_set = ("ER4","LN4")
    end
    for service in service_set
        T = 10000.0
        N = 10000
        record = Record()
        do_experiment(queue, arrival, service, service_rate, T, N, record)
    end
end
