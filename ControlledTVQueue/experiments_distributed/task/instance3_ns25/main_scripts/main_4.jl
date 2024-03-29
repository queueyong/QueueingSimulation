using Distributions, Roots, QuadGK
THIS_FILE_PATH = dirname(@__FILE__)
include("$THIS_FILE_PATH/../../../TVGG1PS_types.jl")
include("$THIS_FILE_PATH/../../../TVGG1PS_functions.jl")

α = 1.0
β = 0.2
γ = 0.001
δ = 2.5e-5

queue = "TVGG1PS" 
control = "SR" 
target = 0.1
arrival = "ER" 
service = "ER" 
fcn = t -> α + β*sin(γ*t)+β*sin(2*γ*t)
fcn_type = "fcn1" 
arrival_table = readdlm("$THIS_FILE_PATH/../../../table/arrival/table_fcn1_arrival.txt")
service_table = readdlm("$THIS_FILE_PATH/../../../table/service/SR/0.1/ER,ER/table_fcn1_service_SR_ER_ER.txt")
T = 22000.0
N = 400
record = Record()
LOG_PATH = "$THIS_FILE_PATH/../log" 

do_experiment_table_slave(queue, control, target, arrival, service, fcn, fcn_type, arrival_table, service_table, T, N, record, LOG_PATH, 4)
