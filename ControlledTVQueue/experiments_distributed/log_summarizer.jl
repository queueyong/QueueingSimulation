mutable struct Instance
    queue::String
    control::String
    target::Float64
    arrival::String
    service::String
    fcn_type::String
    T::Float64
    N::Int64
end

# define instances
instance_set = Instance[]
push!(instance_set, Instance("TVGG1PS","SR",0.1,"EXP","EXP","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",0.1,"EXP","EXP","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","SR",0.1,"ER","ER","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",0.1,"ER","ER","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","SR",0.1,"LN","LN","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",0.1,"LN","LN","fcn1",22000.0,10000))

push!(instance_set, Instance("TVGG1PS","SR",10.0,"EXP","EXP","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",10.0,"EXP","EXP","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","SR",10.0,"ER","ER","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",10.0,"ER","ER","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","SR",10.0,"LN","LN","fcn1",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",10.0,"LN","LN","fcn1",22000.0,10000))

push!(instance_set, Instance("TVGG1PS","SR",0.1,"EXP","EXP","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",0.1,"EXP","EXP","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","SR",0.1,"ER","ER","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",0.1,"ER","ER","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","SR",0.1,"LN","LN","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",0.1,"LN","LN","fcn2",22000.0,10000))

push!(instance_set, Instance("TVGG1PS","SR",10.0,"EXP","EXP","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",10.0,"EXP","EXP","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","SR",10.0,"ER","ER","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",10.0,"ER","ER","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","SR",10.0,"LN","LN","fcn2",22000.0,10000))
push!(instance_set, Instance("TVGG1PS","PD",10.0,"LN","LN","fcn2",22000.0,10000))

function summarize_log(file_stream, LOG_FILE_PATH, LOG_FILE_NAME)
    file_array = readdlm("$LOG_FILE_PATH/$LOG_FILE_NAME")
    nrows, ncols = size(file_array)

    # write time axis
    time_axis = file_array[1,:]
    writedlm(file_stream, transpose(time_axis))

    # calculate avg and ci's
    avg, upper_ci, lower_ci = Float64[], Float64[], Float64[]
    for t in 1:length(time_axis)
        n = nrows-1
        push!(avg, sum(file_array[2:end,t])/n)
        σ = std(file_array[2:end,t])
        push!(upper_ci, avg[t]+1.96*(σ/sqrt(n)))
        push!(lower_ci, avg[t]-1.96*(σ/sqrt(n)))
    end
    writedlm(file_stream, transpose(avg)) # average
    writedlm(file_stream, transpose(upper_ci)) # upper ci
    writedlm(file_stream, transpose(lower_ci)) # lower ci
end

task_set = [(3,25,10),(4,25,7),(9,25,20),(10,25,7),(15,25,10),(16,25,5),(21,25,35),(22,25,10)] # instance: 3,4,9,10,15,16,21,22

THIS_FILE_PATH = dirname(@__FILE__)

for task in task_set
    instance = task[1]
    ns = task[2]
    task_name = "instance$(instance)_ns$(ns)"

    queue = instance_set[instance].queue
    control = instance_set[instance].control
    target = instance_set[instance].target
    arrival = instance_set[instance].arrival
    service = instance_set[instance].service
    fcn_type = instance_set[instance].fcn_type
    T = instance_set[instance].T
    N = instance_set[instance].N

    LOG_FILE_PATH = "$THIS_FILE_PATH/log/$queue/$target/$arrival,$service/$control"
    LOG_FILE_NAME = "$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(N)"
    SUM_LOG_FILE_PATH = "$THIS_FILE_PATH/sum_log/$queue/$target/$arrival,$service/$control"
    SUM_LOG_FILE_NAME = "$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(N)_sum"

    ql_file = open("$SUM_LOG_FILE_PATH/$(SUM_LOG_FILE_NAME)_queue_length.txt", "w")
    st_file = open("$SUM_LOG_FILE_PATH/$(SUM_LOG_FILE_NAME)_sojourn_time.txt", "w")

    summarize_log(ql_file, LOG_FILE_PATH, "$(LOG_FILE_NAME)_queue_length.txt")
    summarize_log(st_file, LOG_FILE_PATH, "$(LOG_FILE_NAME)_sojourn_time.txt")

    close(ql_file)
    close(st_file)
end
