function summarize_log(LOG_FILE::String, SUM_LOG_FILE::String)
    log_file = readdlm(LOG_FILE)


    X = time_axis(file_array)
    avg, upper_ci, lower_ci = Float64[], Float64[], Float64[]
    for i in 1:length(X)
      values = convert(Array{Float64,1}, file_array[2:end,i])
      effective_values = Float64[]
      for v in values
        if v != -1
         push!(effective_values,v)
        end
      end
      n = length(effective_values)
      σ = std(effective_values)
      push!(avg, sum(values)/n)
      push!(upper_ci, avg[i]+1.96*(σ/sqrt(n)))
      push!(lower_ci, avg[i]-1.96*(σ/sqrt(n)))
    end
    writedlm(stream, transpose(X)) # time axis
    writedlm(stream, transpose(avg)) # average
    writedlm(stream, transpose(upper_ci)) # upper ci
    writedlm(stream, transpose(lower_ci)) # lower ci
end


LOG_PATH = "$(dirname(@__FILE__))/log"
SUM_LOG_PATH = "$(dirname(@__FILE__))/sum_log"
queue = "TVGG1PS"
target = 10.0
arrival = "ER"
service = "ER"
fcn_type = "fnc2"
T = 2500.0
N = 10000
control = "SR"

ql_log_file = "$LOG_PATH/$queue/$target/$arrival,$service/$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(N)_queue_length.txt"
st_log_file = "$LOG_PATH/$queue/$target/$arrival,$service/$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(N)_sojourn_time.txt"

ql_log_array = readdlm(ql_log_file)
st_log_array = readdlm(st_log_file)
time_array = st_log_array[1,:]
avg, upper_ci, lower_ci = Float64[], Float64[], Float64[]
for i in 1:length(time_array)
  values = convert(Array{Float64,1}, st_log_array[2:end,i])
  effective_values = Float64[]
  for v in values
#    if v != -1
     push!(effective_values,v)
#    end
  end
  n = length(effective_values)
  σ = std(effective_values)
  push!(avg, sum(values)/n)
  push!(upper_ci, avg[i]+1.96*(σ/sqrt(n)))
  push!(lower_ci, avg[i]-1.96*(σ/sqrt(n)))
end

using PyPlot

plt = PyPlot
plt.figure()
plt.plot(time_array,avg)
plt.ylim(0.0,20.0)


writedlm(stream, transpose(X)) # time axis
writedlm(stream, transpose(avg)) # average
writedlm(stream, transpose(upper_ci)) # upper ci
writedlm(stream, transpose(lower_ci)) # lower ci


for queue in ["TVGG1PS"]
    for target in [0.1]
        for arrival in ["ER"]
            for service in ["ER"]
                for fcn in ["fnc1", "fnc2"]
                    for control in ["SR","PD"]
                        for performance in ["queue_length","sojourn_time"]
                            LOG_FILE = "$LOG_PATH/$queue/$target/$arrival,$service/$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(N)_$(performance).txt"
                            SUM_LOG_FILE = "$SUM_LOG_PATH/$queue/$target/$arrival,$service/$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(fcn_type)_$(T)_$(N)_$(performance).txt"
                            summarize_log(LOG_FILE,SUM_LOG_FILE)
                        end
                    end
                end
            end
        end
    end
end



cd(dirname(Base.source_path()))
queue_length_log_file_path = "../$(queue)/logs/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_queue_length.txt"
sojourn_time_log_file_path = "../$(queue)/logs/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sojourn_time.txt"
queue_length_sum_log_save_path = "../$(queue)/sum_logs/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sum_queue_length.txt"
sojourn_time_sum_log_save_path = "../$(queue)/sum_logs/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sum_sojourn_time.txt"
queue_length_file = open(queue_length_log_file_path)
sojourn_time_file = open(sojourn_time_log_file_path)
queue_length_file_array = readdlm(queue_length_file)
sojourn_time_file_array = readdlm(sojourn_time_file)
close(queue_length_file)
close(sojourn_time_file)
queue_length_sum_log_file = open(queue_length_sum_log_save_path, "w")
sojourn_time_sum_log_file = open(sojourn_time_sum_log_save_path, "w")
save_mean_with_ci(queue_length_file_array, queue_length_sum_log_file)
save_mean_with_ci(sojourn_time_file_array, sojourn_time_sum_log_file)
close(queue_length_sum_log_file)
close(sojourn_time_sum_log_file)
