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

# 합치고자 하는 task 정의
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

    INTEGRATED_FILE_PATH = "$THIS_FILE_PATH/log/$queue/$target/$arrival,$service/$control"
    FILE_NAME = "$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(N)"
    SLAVE_FILE_NAME = "$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(Int(N/ns))"

    ql_file = open("$INTEGRATED_FILE_PATH/$(FILE_NAME)_queue_length.txt","w")
    st_file = open("$INTEGRATED_FILE_PATH/$(FILE_NAME)_sojourn_time.txt","w")

    SLAVE_LOG_PATH = "$THIS_FILE_PATH/task/$task_name/log"
    time = readdlm("$SLAVE_LOG_PATH/$(SLAVE_FILE_NAME)_queue_length_slave1.txt")[1,:]

    # write time axis
    writedlm(ql_file, transpose(time))
    writedlm(st_file, transpose(time))

    # write logs
    for k in 1:ns
        ql_slave_file = "$SLAVE_LOG_PATH/$(SLAVE_FILE_NAME)_queue_length_slave$k.txt"
        st_slave_file = "$SLAVE_LOG_PATH/$(SLAVE_FILE_NAME)_sojourn_time_slave$k.txt"
        ql_slave_file_array = readdlm(ql_slave_file)[2:end,:]
        st_slave_file_array = readdlm(st_slave_file)[2:end,:]
        writedlm(ql_file, ql_slave_file_array)
        writedlm(st_file, st_slave_file_array)
    end
    close(ql_file)
    close(st_file)
end
