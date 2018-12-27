using Distributions, Roots, QuadGK
THIS_FILE_PATH = dirname(@__FILE__)
include("$THIS_FILE_PATH/../../../TVGG1PS_types.jl")
include("$THIS_FILE_PATH/../../../TVGG1PS_functions.jl")

α = 1.0
β = 0.2
γ = 0.001
δ = 2.5e-5

queue = "TVGG1PS" 
control = "PD" 
target = 10.0
arrival = "ER" 
service = "ER" 
fcn = t -> α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t
fcn_type = "fcn2" 
arrival_table = readdlm("$THIS_FILE_PATH/../../../table/arrival/table_fcn2_arrival.txt")
service_table = readdlm("$THIS_FILE_PATH/../../../table/service/PD/10.0/ER,ER/table_fcn2_service_PD_ER_ER.txt")
T = 22000.0
N = 400
record = Record()
LOG_PATH = "$THIS_FILE_PATH/../log" 

do_experiment_table_slave(queue, control, target, arrival, service, fcn, fcn_type, arrival_table, service_table, T, N, record, LOG_PATH, 14)
