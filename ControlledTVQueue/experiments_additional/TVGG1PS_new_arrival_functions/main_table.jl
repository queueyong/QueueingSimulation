using Distributions, Roots, QuadGK

THIS_FILE_PATH = dirname(@__FILE__)
include("$THIS_FILE_PATH/TVGG1PS_types.jl")
include("$THIS_FILE_PATH/TVGG1PS_functions.jl")

T = 22000.0 # 실제 시뮬레이션하는 시간
N = 10

α = 1.0
β = 0.2
γ = 0.001
δ = 0.000025

# fcn
fcn = t -> α + β*sin(γ*t)+β*sin(2*γ*t)
fcn = t -> α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t
fcn_type = "fcn2"

queue = "TVGG1PS"
control = "SR"
target = 10.0
arrival = "ER"
service = "ER"

record = Record()
LOG_PATH = "$THIS_FILE_PATH/log"
arrival_table = readdlm("$THIS_FILE_PATH/table/arrival/table_$(fcn_type)_arrival.txt")
service_table = readdlm("$THIS_FILE_PATH/table/service/$control/$target/$arrival,$service/table_$(fcn_type)_service_$(control)_$(arrival)_$(service).txt")
do_experiment_table(queue, control, target, arrival, service, fcn, fcn_type, arrival_table, service_table, T, N, record, LOG_PATH)
