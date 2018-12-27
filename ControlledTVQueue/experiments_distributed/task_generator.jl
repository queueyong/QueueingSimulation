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

# define functions
α = 1.0
β = 0.2
γ = 0.001
δ = 0.000025

fcn1 = t -> α + β*sin(γ*t)+β*sin(2*γ*t)
fcn1_string = "t -> α + β*sin(γ*t)+β*sin(2*γ*t)"
fcn2 = t -> α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t
fcn2_string = "t -> α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t"
fcn_dict = Dict("fcn1"=>fcn1_string, "fcn2"=>fcn2_string)

## 요렇게 다 만들고싶다
# (instance번호, slave수, 1rep 시간)
task_set = [(3,25,10),(4,25,7),(9,25,20),(10,25,7),(15,25,10),(16,25,5),(21,25,35),(22,25,10)] # instance: 3,4,9,10,15,16,21,22
# single_rep_expected_time_set = [10,7,20,7,10,5,35,10]
#######################################################################################################################

## generate all
for pair in task_set
    instance = pair[1]
    num_slave = pair[2]
    task_name = "instance$(instance)_ns$(num_slave)"
    THIS_FILE_PATH = dirname(@__FILE__)

    queue = instance_set[instance].queue
    control = instance_set[instance].control
    target = instance_set[instance].target
    arrival = instance_set[instance].arrival
    service = instance_set[instance].service
    fcn_type = instance_set[instance].fcn_type
    T = instance_set[instance].T
    N = instance_set[instance].N

    bebop_time_hour = ceil(pair[3]*(N/num_slave)/3600)

    TASK_PATH = "$THIS_FILE_PATH/task"
    cd(TASK_PATH)
    mkdir("$task_name")
    mkdir("$task_name/log")
    mkdir("$task_name/main_scripts")
    mkdir("$task_name/shell_scripts")

    #### script begins
    for k in 1:num_slave
        # generate main_scripts
        file = open("$TASK_PATH/$task_name/main_scripts/main_$k.jl","w")
        println(file, "using Distributions, Roots, QuadGK")
        println(file, "THIS_FILE_PATH = dirname(@__FILE__)")
        println(file, """include("\$THIS_FILE_PATH/../../../TVGG1PS_types.jl")""")
        println(file, """include("\$THIS_FILE_PATH/../../../TVGG1PS_functions.jl")""")
        println(file, "")
        println(file, "α = $α")
        println(file, "β = $β")
        println(file, "γ = $γ")
        println(file, "δ = $δ")
        println(file, "")
        println(file, """queue = "$queue" """)
        println(file, """control = "$control" """)
        println(file, "target = $target")
        println(file, """arrival = "$arrival" """)
        println(file, """service = "$service" """)
        println(file, "fcn = $(fcn_dict[fcn_type])")
        println(file, """fcn_type = "$(fcn_type)" """)
        println(file, """arrival_table = readdlm("\$THIS_FILE_PATH/../../../table/arrival/table_$(fcn_type)_arrival.txt")""")
        println(file, """service_table = readdlm("\$THIS_FILE_PATH/../../../table/service/$control/$target/$arrival,$service/table_$(fcn_type)_service_$(control)_$(arrival)_$(service).txt")""")
        println(file, "T = $T")
        println(file, "N = $(Int(N/num_slave))")
        println(file, "record = Record()")
        println(file, """LOG_PATH = "\$THIS_FILE_PATH/../log" """)
        println(file, "")
        println(file, "do_experiment_table_slave(queue, control, target, arrival, service, fcn, fcn_type, arrival_table, service_table, T, N, record, LOG_PATH, $k)")
        close(file)

        # generate shell_scripts
        file = open("$TASK_PATH/$task_name/shell_scripts/shell_$k.jl","w")
        println(file, "\#!/bin/bash")
        println(file, "\#SBATCH -J SL$(instance)ns$(k)")
        println(file, "\#SBATCH -p bdwall")
        println(file, "\#SBATCH -A NEXTGENOPT")
        println(file, "\#SBATCH -N 1")
        println(file, "\#SBATCH -t $(Int(bebop_time_hour)):00:00")
        println(file, "")
        println(file, "julia $TASK_PATH/$task_name/main_scripts/main_$k.jl")
        close(file)
    end
end


#=
#### script begins
k = 1

file = open("$TASK_PATH/$task_name/main_scripts/main_$k.jl","w")
println(file, "using Distributions, Roots, QuadGK")
println(file, "THIS_FILE_PATH = dirname(@__FILE__)")
println(file, """include("\$THIS_FILE_PATH/../../../TVGG1PS_types.jl")""")
println(file, """include("\$THIS_FILE_PATH/../../../TVGG1PS_functions.jl")""")
println(file, "")
println(file, "α = $α")
println(file, "β = $β")
println(file, "γ = $γ")
println(file, "δ = $δ")
println(file, "")
println(file, """queue = "$queue" """)
println(file, """control = "$control" """)
println(file, "target = $target")
println(file, """arrival = "$arrival" """)
println(file, """service = "$service" """)
println(file, "fcn = $(fcn_dict[fcn_type])")
println(file, """fcn_type = "$(fcn_type)" """)
println(file, """arrival_table = readdlm("\$THIS_FILE_PATH/../../../table/arrival/table_$(fcn_type)_arrival.txt")""")
println(file, """service_table = readdlm("\$THIS_FILE_PATH/../../../table/service/$control/$target/$arrival,$service/table_$(fcn_type)_service_$(control)_$(arrival)_$(service).txt")""")
println(file, "T = $T")
println(file, "N = $(Int(N/num_slave))")
println(file, "record = Record()")
println(file, """LOG_PATH = "\$THIS_FILE_PATH/../log" """)
println(file, "")
println(file, "do_experiment_table_slave(queue, control, target, arrival, service, fcn, fcn_type, arrival_table, service_table, T, N, record, LOG_PATH, $k)")
close(file)

=#
