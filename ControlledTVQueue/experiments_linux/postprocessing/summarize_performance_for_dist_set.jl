using DataFrames, CSV

function getAmplitudePerformance(file_array::Matrix{Float64}, period::Float64)
    T = file_array[1,end]
    i = size(file_array)[2]
    while file_array[1,i] >= T - period
        i -= 1
    end
    max = maximum(file_array[2,i:end])
    min = minimum(file_array[2,i:end])
    amplitude = (max - min)/2
    relative_amplitude = amplitude/getNumericalAverage(file_array, period)*100
    return amplitude , relative_amplitude
end

function getTargetPerformance(file_array::Matrix{Float64}, period::Float64, real_target::Float64)
    controlled_target = getNumericalAverage(file_array, period)
    target_gap = (controlled_target - real_target)/real_target
    return controlled_target, target_gap
end

function getNumericalAverage(file_array::Matrix{Float64}, period::Float64)

    interval = file_array[1,2]
    T = file_array[1,end]
    area = 0.0
    i = size(file_array)[2]
    while file_array[1,i] >= T - period
        area += file_array[2,i]*interval
        i -= 1
    end

    return area/period
end


# read sum_log file

# calculate performance measure

# write on a file (.cvs)

N = 10000
queue_set = ("TVGG1", "TVGG1PS")
dist_set = (("EXP","EXP") , ("LN","LN") , ("ER","ER") , ("LN","ER") , ("ER","LN"))
target_set = (0.1, 1.0, 10.0)
gamma_set = (0.001, 0.01, 0.1)
T_set = (20000.0, 2000.0, 2000.0)
control_set = ("SR", "PD")

for queue in queue_set
    for d in dist_set
        arrival, service = d
        df = DataFrame(target=Any[], gamma=Any[], SR_amplitude=Any[], SR_average=Any[], PD_amplitude=Any[], PD_average=Any[])
        for s in target_set
            for (γ,T) in zip(gamma_set, T_set)
                period = 2*π/γ
                tempv = Any[]
                push!(tempv, s)
                push!(tempv, γ)
                for control in control_set
                    sojourn_time_file_path = "$(dirname(@__FILE__))/../sum_logs/$(queue)/$(s)/$(arrival),$(service)/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(γ)_$(T)_$(N)_sum_sojourn_time.txt"
                    sojourn_time_file = open(sojourn_time_file_path)
                    sojourn_time_file_array = readdlm(sojourn_time_file)
                    close(sojourn_time_file)
                    push!(tempv, "$(round(getAmplitudePerformance(sojourn_time_file_array,period)[1],4)) ($(round(getAmplitudePerformance(sojourn_time_file_array,period)[2],2))%)")
                    push!(tempv, "$(round(getTargetPerformance(sojourn_time_file_array,period,s)[1],4)) ($(round(getTargetPerformance(sojourn_time_file_array,period,s)[2],2))%)")
                end
                push!(df,tempv)
            end
        end
        CSV.write("$(dirname(@__FILE__))/../performance_analysis/$(queue)_$(arrival)_$(service)_performance_analysis.csv" , df)
    end
end

T = 20000.0
N = 10000
queue = "TVGG1"
d = dist_set[1]
arrival, service = d
s = target_set[1]
γ = gamma_set[1]
control = control_set[1]


cd(dirname(Base.source_path()))
pwd()
queue_length_file_path = "../sum_logs/$(queue)/$(s)/$(arrival),$(service)/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(γ)_$(T)_$(N)_sum_queue_length.txt"
sojourn_time_file_path = "../sum_logs/$(queue)/$(s)/$(arrival),$(service)/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(γ)_$(T)_$(N)_sum_sojourn_time.txt"
queue_length_file = open(queue_length_file_path)
sojourn_time_file = open(sojourn_time_file_path)
queue_length_file_array = readdlm(queue_length_file)
sojourn_time_file_array = readdlm(sojourn_time_file)
close(queue_length_file)
close(sojourn_time_file)

interval = queue_length_file_array[1,2]
T = queue_length_file_array[1,end]
period = 2*π/γ
area = 0.0
i = size(queue_length_file_array)[2]
while queue_length_file_array[1,i] >= T - period
    area += queue_length_file_array[2,i]*interval
    i -= 1
end

area/period
