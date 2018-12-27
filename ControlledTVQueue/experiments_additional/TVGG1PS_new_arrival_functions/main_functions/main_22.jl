using Distributions, Roots, QuadGK

THIS_FILE_PATH = dirname(@__FILE__)
include("$THIS_FILE_PATH/../TVGG1PS_types.jl")
include("$THIS_FILE_PATH/../TVGG1PS_functions.jl")

T = 21000.0 # 실제 시뮬레이션하는 시간
N = 10000

α = 1.0
β = 0.2
γ = 0.001
δ = 0.000025

queue = "TVGG1PS"
control = "PD"
target = 10.0
arrival = "ER"
service = "ER"
fcn = t -> α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t
period = T
fcn_type = "fcn2"
record = Record()
LOG_PATH = "$THIS_FILE_PATH/../log/$queue/$target/$arrival,$service"
do_experiment(queue, control, target, arrival, service, fcn, period, fcn_type, T, N, record, LOG_PATH)
