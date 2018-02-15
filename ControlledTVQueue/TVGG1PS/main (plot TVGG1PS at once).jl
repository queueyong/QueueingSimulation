cd(dirname(Base.source_path()))
using PyPlot

function time_axis(file_array::Matrix{Any})
 temp_array = file_array[1,:]
 x = Float64[]
 i = 1
 while typeof(temp_array[i]) != SubString{String}
  push!(x, temp_array[i])
  i += 1
 end
 return x
end

function time_axis(file_array::Matrix{Float64})
 temp_array = file_array[1,:]
 return temp_array
end

function plot_tail_probatility(file_array::Any, plt::Module)
  δ_set = file_array[1,:]
  x = file_array[2,:]
  for i in 3:(size(file_array)[1]-1)
    y = file_array[i,:]
    plt.plot(x,y,label="P[W(t)>$(δ_set[i-2])]")
  end
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

function plot_mean_number(file_array::Any, plt::Module, _label::String, _color::String)
  x = file_array[1,:]
  y = file_array[2,:]
  plt.plot(x,y,color=_color, label=_label)
end

function plot_virtual_sojourn_time(file_array::Any, plt::Module, _label::String, _color::String)
  x = file_array[1,:]
  y = file_array[2,:]
  plt.plot(x,y,color=_color, label=_label)
end

function plot_arrival_rate_function(x::Any, coeff::Tuple, plt::Module)
 λ = t -> coeff[1] + coeff[2]*sin(coeff[3]*t)
 y = λ(x)
 plt.plot(x,y,color="darkblue",linestyle="--", label="arrival rate")
end


# ALL AT ONCE
queue = "TVGG1PS"
queue_type = L"$G_t/G_t/1/PS$"

control_set = ("RM", "SR", "PSA", "DM")
param_letter_set = ("ρ", "ξ", "ζ", "δ")
param_set = (0.8,0.2,1.0, 0.25)

distribution_set = [("Exponential","Exponential"), ("Erlang","Erlang"), ("Hyperexponential","Hyperexponential"), ("Hyperexponential","Erlang")]
gamma_set = (0.001, 0.01, 0.1, 1.0 , 10.0)
time_set = (20000.0, 2000.0, 2000.0, 2000.0, 2000.0)
colors = ("lime","tomato","darkblue","purple")

for c in 1:length(control_set)
 for d in 1:length(distribution_set)
   plt = PyPlot
   plt.figure(figsize=(10.5,15))
   arrival = distribution_set[d][1]
   service = distribution_set[d][2]
   param_letter = param_letter_set[c]
   param = param_set[c]
   control = control_set[c]
   for i in 1:length(gamma_set)
     plt.subplot(5,1,i)
     coeff = (1.0, 0.2, gamma_set[i])
     time = time_set[i]
     plt.title("$queue_type, $arrival/$service, λ(t) = $(coeff[1])+$(coeff[2])sin($(coeff[3])t), $control ($param_letter = $param)")
     plt.xlabel("time")
     plt.xlim(0,time)
     num_in_queue_file = open("C:/Users/YOC/Google 드라이브/YOC-Data/Summarized Simulation Log/$queue/$control/$arrival, $service/(sum) num in queue ($queue, gamma $(coeff[3]), control $control, arrival $arrival, service $service, time $time).txt")
     virtual_sojourn_time_file = open("C:/Users/YOC/Google 드라이브/YOC-Data/Summarized Simulation Log/$queue/$control/$arrival, $service/(sum) virtual sojourn time ($queue, gamma $(coeff[3]), control $control, arrival $arrival, service $service, time $time).txt")
     num_in_queue_file_array = readdlm(num_in_queue_file) # convert .txt file to array
     virtual_sojourn_time_file_array = readdlm(virtual_sojourn_time_file)
     close(num_in_queue_file)
     close(virtual_sojourn_time_file)
     x = linspace(0.0, time, 10000)
     plot_mean_number(num_in_queue_file_array, plt)
     plot_ci(num_in_queue_file_array, plt)
     plot_virtual_sojourn_time(virtual_sojourn_time_file_array, plt)
     plot_ci(virtual_sojourn_time_file_array, plt)
     plot_arrival_rate_function(x, coeff, plt)
     plt.legend(loc = "right", fontsize = 8)
   end
   plt.tight_layout()
   plt.savefig("./plot type 1 (mean number and sojourn time) ($queue, $control, $arrival, $service).pdf")
   close()

   plt.figure(figsize=(10.5,15))
   for i in 1:length(gamma_set)
     plt.subplot(5,1,i)
     coeff = (1.0, 0.2, gamma_set[i])
     time = time_set[i]
     plt.title("$queue_type, $arrival/$service, λ(t) = $(coeff[1])+$(coeff[2])sin($(coeff[3])t), $control ($param_letter = $param)")
     plt.xlabel("time")
     plt.xlim(0,time)
     tail_pr_file = open("C:/Users/YOC/Google 드라이브/YOC-Data/Summarized Simulation Log/$queue/$control/$arrival, $service/(sum) tail probability ($queue, gamma $(coeff[3]), control $control, arrival $arrival, service $service, time $time).txt")
     tail_pr_file_array = readdlm(tail_pr_file)
     close(tail_pr_file)
     x = linspace(0.0, time, 10000)
     plot_tail_probatility(tail_pr_file_array, plt)
     plot_arrival_rate_function(x, coeff, plt)
     plt.legend(loc = "right", fontsize = 8)
   end
   plt.tight_layout()
   plt.savefig("./plot type 2 (tail probability) ($queue, $control, $arrival, $service).pdf")
   close()
 end

 plt = PyPlot
 plt.figure(figsize=(10.5,15))
 param_letter = param_letter_set[c]
 param = param_set[c]
 control = control_set[c]
 for i in 1:length(gamma_set)
   γ = gamma_set[i]
   time = time_set[i]
   x = linspace(0.0, time, 10000)
   ## plotting
   plt.subplot(5,1,i)
   coeff = (1.0, 0.2, γ)
   plt.title("$queue_type, λ(t) = $(coeff[1])+$(coeff[2])sin($(coeff[3])t), $control ($param_letter = $param)")
   plt.xlabel("time")
   plt.xlim(0,time)
   labels = ("E[Q(t)] (EXP/EXP)","E[Q(t)] (H2/H2)","E[Q(t)] (H2/E2)","E[Q(t)] (E2/E2)")
   for j in 1:length(distribution_set)
     num_in_queue_file = open("C:/Users/YOC/Google 드라이브/YOC-Data/Summarized Simulation Log/$queue/$control/$(distribution_set[j][1]), $(distribution_set[j][2])/(sum) num in queue ($queue, gamma $(coeff[3]), control $control, arrival $(distribution_set[j][1]), service $(distribution_set[j][2]), time $time).txt")
     num_in_queue_file_array = readdlm(num_in_queue_file) # convert .txt file to array
     close(num_in_queue_file)
     plot_mean_number(num_in_queue_file_array, plt, labels[j], colors[j])
     plot_ci(num_in_queue_file_array, plt)
   end
   plot_arrival_rate_function(x, coeff, plt)
   plt.legend(loc = "upper right", fontsize = 8)
 end
 plt.tight_layout()
 plt.savefig("./plot type 3 (mean number) ($queue, $control, general distributions).pdf")
 close()

 plt.figure(figsize=(10.5,15))
 param_letter = param_letter_set[c]
 param = param_set[c]
 control = control_set[c]
 for i in 1:length(gamma_set)
   γ = gamma_set[i]
   time = time_set[i]
   x = linspace(0.0, time, 10000)
   ## plotting
   plt.subplot(5,1,i)
   coeff = (1.0, 0.2, γ)
   plt.title("$queue_type, λ(t) = $(coeff[1])+$(coeff[2])sin($(coeff[3])t), $control ($param_letter = $param)")
   plt.xlabel("time")
   plt.xlim(0,time)
   labels = ("E[S(t)] (EXP/EXP)","E[S(t)] (H2/H2)","E[S(t)] (H2/E2)","E[S(t)] (E2/E2)")
   for j in 1:length(distribution_set)
     virtual_sojourn_time_file = open("C:/Users/YOC/Google 드라이브/YOC-Data/Summarized Simulation Log/$queue/$control/$(distribution_set[j][1]), $(distribution_set[j][2])/(sum) virtual sojourn time ($queue, gamma $(coeff[3]), control $control, arrival $(distribution_set[j][1]), service $(distribution_set[j][2]), time $time).txt")
     virtual_sojourn_time_file_array = readdlm(virtual_sojourn_time_file)
     close(virtual_sojourn_time_file)
     plot_virtual_sojourn_time(virtual_sojourn_time_file_array, plt, labels[j], colors[j])
     plot_ci(virtual_sojourn_time_file_array, plt)
   end
   plot_arrival_rate_function(x, coeff, plt)
   plt.legend(loc = "upper right", fontsize = 8)
 end
 plt.tight_layout()
 plt.savefig("./plot type 4 (mean sojourn) ($queue, $control, general distributions).pdf")
 close()
end
