using Distributions, Roots, QuadGK

THIS_FILE_PATH = dirname(@__FILE__)
include("$THIS_FILE_PATH/TVGG1PS_types.jl")
include("$THIS_FILE_PATH/TVGG1PS_functions.jl")

T = 20000.0
N = 10
fcn2(t) = 0.00004*t*sin(6*π/T*t)+0.00004*t*sin(2*6*π/T*t)+2*0.00004*t+0.1
record = Record()
do_experiment("TVGG1PS", "SR", 0.1, "ER", "ER", fcn2, T/3, "fcn2", T, N, record, "$THIS_FILE_PATH/log")


# main
queue = "TVGG1PS"
N = 10000
control_set = ["PD", "SR"]
#dist_set = [ ("EXP","EXP"), ("LN","LN"), ("ER","ER"), ("ER","LN") , ("LN","ER")]
dist_set = [ ("ER2","ER2"), ("LN2","LN2") ]
target_set = [0.1]
γ_set = [0.1, 0.01, 0.001]
T_set = [2000.0, 2000.0, 20000.0]
γT_zip = zip(γ_set, T_set)


for control in control_set
    for dist in dist_set
        for target in target_set
            for γT in γT_zip
                arrival = dist[1]
                service = dist[2]
                coeff = (1.0, 0.2, γT[1])
                T = γT[2]   # run time
                record = Record()
                do_experiment(queue, control, target, arrival, service, coeff, T, N, record, LOG_PATH)
            end
        end
    end
end
println("Done")

#=
for control in control_set
    for dist in dist_set
        for target in target_set
            for i in 1:3
                arrival = dist[1]
                service = dist[2]
                coeff = (1.0, 0.2, γ_set[i])
                T = T_set[i]
                record = Record()
                do_experiment(queue, control, target, arrival, service, coeff, T, N, record)
            end
        end
    end
end
=#

#=
queue = "TVGG1PS"
N = 10000
arrival = "ER"
service = "ER"
coeff = (1.0, 0.2, 0.1)
T = 2000.0
record = Record()
do_experiment(queue, control, target, arrival, service, coeff, T, N, record)
=#
