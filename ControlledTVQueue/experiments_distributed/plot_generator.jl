using PyPlot

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
fcn_function_dict = Dict("fcn1"=>fcn1, "fcn2"=>fcn2)

function getNumericalAverage(array)
    sum = 0.0
    for a in array
        sum += a
    end
    return sum/length(array)
end

function getFunctionValues(x, fcn)
    temp = Float64[]
    for val in x
        push!(temp, fcn(val))
    end
    return temp
end

task_set = [(3,25,10),(4,25,7),(9,25,20),(10,25,7),(15,25,10),(16,25,5),(21,25,35),(22,25,10)] # instance: 3,4,9,10,15,16,21,22

THIS_FILE_PATH = dirname(@__FILE__)

for task in task_set
    instance = task[1]

    queue = instance_set[instance].queue
    control = instance_set[instance].control
    target = instance_set[instance].target
    arrival = instance_set[instance].arrival
    service = instance_set[instance].service
    fcn_type = instance_set[instance].fcn_type
    T = instance_set[instance].T
    N = instance_set[instance].N

    PLOT_SAVE_DIR = "$THIS_FILE_PATH/plot"
    SUM_LOG_FILE_PATH = "$THIS_FILE_PATH/sum_log/$queue/$target/$arrival,$service/$control"
    SUM_LOG_FILE_NAME = "$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(N)_sum"

    ql_sum_array = readdlm("$SUM_LOG_FILE_PATH/$(SUM_LOG_FILE_NAME)_queue_length.txt")
    st_sum_array = readdlm("$SUM_LOG_FILE_PATH/$(SUM_LOG_FILE_NAME)_sojourn_time.txt")

    time = ql_sum_array[1,:]
    #fig = figure(figsize=(10,5))
    fig = figure()
    ax1 = fig[:add_subplot](1,1,1)
    ax1[:locator_params](nbins=4)
    ax1[:set_xlim](0.0,20000.0)

    if target == 0.1
        ax1[:set_ylim](0.0,3*target)
        #ax1[:set_ylim](0.0,3*getNumericalAverage(st_sum_array[2,:]))
    else
        ax1[:set_ylim](0.0,2.5*target)
        #ax1[:set_ylim](0.0,2*getNumericalAverage(st_sum_array[2,:]))
    end

    fs = 12

    setp(ax1[:get_xticklabels](),fontsize=fs)
    setp(ax1[:get_yticklabels](),fontsize=fs)
    ax1[:plot](time,ql_sum_array[2,:], label="E[Q(t)]", color="lime")
    ax1[:plot](time,ql_sum_array[3,:], color="black", linestyle = "--", linewidth = 0.1)
    ax1[:plot](time,ql_sum_array[4,:], color="black", linestyle = "--", linewidth = 0.1)
    ax1[:plot](time,st_sum_array[2,:], label="E[R(t)]", color="tomato")
    ax1[:plot](time,st_sum_array[3,:], color="black", linestyle = "--", linewidth = 0.1)
    ax1[:plot](time,st_sum_array[4,:], color="black", linestyle = "--", linewidth = 0.1)
    ax1[:plot](time, -100*time, label="arrival rate", linestyle = "--", color="darkblue")
    ax1[:legend](loc = "upper left", fontsize=fs)
    xlabel("time", fontsize=fs)

    ax2 = ax1[:twinx]()
    ax2[:locator_params](nbins=6)
    ax2[:set_ylim](0.0,5)
    ax2[:plot](time,getFunctionValues(time,fcn_function_dict[fcn_type]), color="darkblue", linestyle = "--")
    setp(ax2[:get_yticklabels](),color="darkblue", fontsize=fs)
    ylabel("λ(t)", color="darkblue", fontsize=fs)
    tight_layout()

    if target == 0.1
        PLOT_FILE_NAME = "$(queue)_01_$(control)_$(arrival)_$(service)_$(fcn_type)_$(Int(T))_$(N)"
    else
        PLOT_FILE_NAME = "$(queue)_10_$(control)_$(arrival)_$(service)_$(fcn_type)_$(Int(T))_$(N)"
    end
    savefig("$PLOT_SAVE_DIR/$PLOT_FILE_NAME.pdf")

end
