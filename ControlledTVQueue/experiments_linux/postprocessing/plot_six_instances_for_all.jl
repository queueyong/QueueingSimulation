cd(dirname(Base.source_path()))
using PyPlot

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
  plt.plot(x,upper_ci,color="black", linestyle="--", linewidth = 0.1)
  plt.plot(x,lower_ci,color="black", linestyle="--", linewidth = 0.1)
end

function plot_arrival_rate_function(x::Any, coeff::Tuple, plt::Module)
 λ = t -> coeff[1] + coeff[2]*sin(coeff[3]*t)
 y = λ(x)
 plt.plot(x,y,color="darkblue",linestyle="--", label="arrival rate")
end

function full_name(short_name::String)
 if short_name == "ER"
  return "Erlang"
 elseif short_name == "LN"
  return "Lognormal"
 elseif short_name == "TVGG1"
  return L"$G_t/G_t/1$"
 elseif short_name == "TVGG1PS"
  return L"$G_t/G_t/1/PS$"
 elseif short_name == "SR"
  return "Square-root control"
 elseif short_name == "PD"
  return "Difference-matching control"
 elseif short_name == "EXP"
  return "Exp"
 end
end


## main part
cd(dirname(Base.source_path()))

queue_set = ("TVGG1", "TVGG1PS")
s_set = (0.1, 1.0, 10.0)
dist_set = (("ER","LN"), ("LN","ER"), ("EXP","EXP"), ("ER","ER"), ("LN","LN"))
N = 10000
control_set = ("SR","PD")
coeff_set = ((1.0, 0.2, 0.001) , (1.0, 0.2, 0.01) , (1.0, 0.2, 0.1))

ylim_factor = 0.0

for queue in queue_set
 for s in s_set
  for dist in dist_set
   arrival, service = dist
   plt = PyPlot
   plt.figure(figsize=(8,10))
   suptitle("$(full_name(queue)), Target sojourn time = $s, Arrival=$(full_name(arrival)), Service=$(full_name(service))")
   for j in 1:6
    plt.subplot(3,2,j)
    if j == 1 || j == 2
     plt.title("$(full_name(control_set[j]))")
    end
    plt.xlabel("time")

    k = 0
    if 1 <= j <= 2
     k = 1
    elseif 3 <= j <=4
     k = 2
    elseif 5 <= j <= 6
     k = 3
    end

    coeff = coeff_set[k]
    control = control_set[isodd(j) ? 1:2]

    T = 0.0
    if j <= 2
     T = 20000.0
    else
     T = 2000.0
    end

    queue_length_file_path = "../sum_logs/$queue/$s/$arrival,$service/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sum_queue_length.txt"
    sojourn_time_file_path = "../sum_logs/$queue/$s/$arrival,$service/$(queue)_$(s)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sum_sojourn_time.txt"
    queue_length_file = open(queue_length_file_path)
    sojourn_time_file = open(sojourn_time_file_path)
    queue_length_file_array = readdlm(queue_length_file)
    sojourn_time_file_array = readdlm(sojourn_time_file)
    close(queue_length_file)
    close(sojourn_time_file)

    if j == 1 || j == 3 || j ==5
     plt.ylabel("γ=$(coeff[3])")
     ylim_factor = getNumericalAverage(sojourn_time_file_array, 2*π/coeff[3])
    end

    plt.xlim(0,T)
    plt.ylim(0,ylim_factor*2)

    X = linspace(0.0,T,10000)
    plot_mean_number(queue_length_file_array, plt)
    plot_ci(queue_length_file_array, plt)
    plot_virtual_sojourn_time(sojourn_time_file_array, plt)
    plot_ci(sojourn_time_file_array, plt)
    #plot_arrival_rate_function(X, coeff, plt)
    plt.legend(loc = "upper right", fontsize = 8)
    #plt.tight_layout()
    plt.savefig("../plots/$(queue)_$(s)_$(arrival)_$(service).pdf")
   end
   plt.close()
  end
 end
end
