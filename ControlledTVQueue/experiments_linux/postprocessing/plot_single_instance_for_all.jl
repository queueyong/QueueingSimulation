using PyPlot

function time_axis(file_array::Any)
 temp_array = file_array[1,:]
 X = Float64[]
 for i in 1:1001
  push!(X, temp_array[i])
 end
 return X
end

function plot_mean_number(file_array::Any, plt::Module)
  x = file_array[1,:]
  y = file_array[2,:]
  plt.plot(x,y,color="lime", label="E[Q(t)]")
end

function plot_virtual_sojourn_time(file_array::Any, plt::Module)
  x = file_array[1,:]
  y = file_array[2,:]
  plt.plot(x,y,color="tomato", label="E[S(t)]")
end

function plot_ci(file_array::Any, plt::Module)
  x = file_array[1,:]
  upper_ci = file_array[3,:]
  lower_ci = file_array[4,:]
  plt.plot(x,upper_ci,color="black", linestyle="--", linewidth = 0.3)
  plt.plot(x,lower_ci,color="black", linestyle="--", linewidth = 0.3)
end

function plot_arrival_rate_function(x::Any, coeff::Tuple, plt::Module)
 λ = t -> coeff[1] + coeff[2]*sin(coeff[3]*t)
 y = λ(x)
 plt.plot(x,y,color="darkblue",linestyle="--", label="arrival rate")
end


cd(dirname(Base.source_path()))
control_set = ("PD",)
coeff_set = ((1.0, 0.2, 0.1) , (1.0, 0.2, 0.01) , (1.0, 0.2, 0.001))
queue_set = ("TVGG1PS",)
s_set = (0.1, 1.0)
dist_set = (("EXP","EXP"), ("ER","ER"), ("LN","LN"), ("ER","LN"), ("LN","ER"))
T_set = (2000.0, 2000.0, 20000.0)
N = 10000

for queue in queue_set
 for s in s_set
  for dist in dist_set
   for control in control_set
    for (coeff, T) in zip(coeff_set,T_set)
        arrival,service = dist[1], dist[2]
        plt = PyPlot
        plt.figure()
        plt.title("$queue, $arrival/$service, λ(t) = $(coeff[1])+$(coeff[2])sin($(coeff[3])t), $control, s = $s")
        plt.xlabel("time")
        plt.xlim(0,T)
        queue_length_file_path = "../sum_logs/$queue/$s/$arrival,$service/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sum_queue_length.txt"
        sojourn_time_file_path = "../sum_logs/$queue/$s/$arrival,$service/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sum_sojourn_time.txt"
        queue_length_file = open(queue_length_file_path)
        sojourn_time_file = open(sojourn_time_file_path)
        queue_length_file_array = readdlm(queue_length_file)
        sojourn_time_file_array = readdlm(sojourn_time_file)
        close(queue_length_file)
        close(sojourn_time_file)

        X = linspace(0.0,T,10000)
        plot_mean_number(queue_length_file_array, plt)
        plot_ci(queue_length_file_array, plt)
        plot_virtual_sojourn_time(sojourn_time_file_array, plt)
        plot_ci(sojourn_time_file_array, plt)
        #plot_arrival_rate_function(X, coeff, plt)
        plt.legend(loc = "right", fontsize = 8)
        plt.savefig("../plots/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3]).pdf")
        plt.close()

    end
   end
  end
 end
end
