cd(dirname(Base.source_path()))
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
  plt.plot(x,y,color="tomato", label="E[S(t)] simulation")
end

function plot_ci(file_array::Any, plt::Module)
  x = file_array[1,:]
  upper_ci = file_array[3,:]
  lower_ci = file_array[4,:]
  plt.plot(x,upper_ci,color="black", linestyle="--", linewidth = 0.1)
  plt.plot(x,lower_ci,color="black", linestyle="--", linewidth = 0.1)
end

function plot_sojourn_time_approximation(X::Any, arrival::String, service::String, μ::Float64, plt::Module)
 #x = linspace(0.0,10000.0,10000)
 λ = 1/mean(arrival)
 β = mean(service)
 ρ = λ*β/μ
 approx = ((scv(arrival)+scv(service))/(1+scv(service)))*(1/(1-ρ))*(β/μ)
 S = [approx for i in 1:length(X)]
 plt.plot(X,S,color="green",linestyle="--", label=L"$E[S(t)]\approx V_{PS}\cdot\frac{1}{1-\rho}\cdot\frac{\beta}{\mu}$")
end

function plot_sojourn_time_approximation_ρ(X::Any, arrival::String, service::String, μ::Float64, plt::Module)
 #x = linspace(0.0,10000.0,10000)
 λ = 1/mean(arrival)
 β = mean(service)
 ρ = λ*β/μ
 approx = ((scv(arrival)+scv(service))/(1+scv(service)))*(ρ/(1-ρ))*(β/μ)
 S = [approx for i in 1:length(X)]
 plt.plot(X,S,color="darkblue",linestyle="--", label=L"$E[S(t)]\approx V_{PS}\cdot\frac{\rho}{1-\rho}\cdot\frac{\beta}{\mu}$")
end

function plot_sojourn_time_approximation_GG1(X::Any, arrival::String, service::String, μ::Float64, plt::Module)
 #x = linspace(0.0,10000.0,10000)
 λ = 1/mean(arrival)
 β = mean(service)
 ρ = λ*β/μ
 approx = (β/μ)+(β/μ)*(ρ/(1-ρ))*((scv(arrival)+scv(service))/2)
 S = [approx for i in 1:length(X)]
 plt.plot(X,S,color="brown",linestyle="--", label=L"$E[S(t)]\approx V\cdot\frac{\rho}{1-\rho}\cdot\frac{\beta}{\mu}+\frac{\beta}{\mu}$")
end

#=
function plot_queue_length_approximation(arrival::String, service::String, μ::Float64, plt::Module)
 x = linspace(0.0,10000.0,10000)
 Q = [ for i in 1:length(x)]
end
=#
function scv(short_dist::String)
 if short_dist ∈ ("ER0","ER1","ER2")
   return 0.5
 elseif short_dist ∈ ("LN0","LN1","LN2")
   return 2.0
 end
end

function Base.mean(short_dist::String)
 if short_dist ∈ ("ER1","LN1")
  return 0.9
 elseif short_dist ∈ ("ER2","LN2")
  return 1.8
 elseif short_dist ∈ ("ER0","LN0")
  return 1.0
 end
end

function full_name(short_name::String)
 if short_name == "ER0"
  return "Erlang (mean=1)"
 elseif short_name == "LN0"
  return "Lognormal (mean=1)"
 elseif short_name == "ER1"
  return "Erlang (mean=0.9)"
 elseif short_name == "LN1"
  return "Lognormal (mean=0.9)"
 elseif short_name == "ER2"
  return "Erlang (mean=1.8)"
 elseif short_name == "LN2"
  return "Lognormal (mean=1.8)"
 elseif short_name == "GG1PS"
  return L"$G/G/1/PS$"
 end
end

# modify here =====================
queue = "GG1PS"
μ_set = (1.0, 2.0)
arrival_set = ("ER0","LN0")
service_set = (("ER1","LN1"),("ER2","LN2"))

#ylim_factor = 1.5
#======================#
T = 10000.0
N = 10000
ρ = 0.9
plt = PyPlot
plt.figure(figsize=(16,11))

suptitle("$(full_name(queue)), ρ = λβ/μ = 0.9, SCV: Erlang ⇒ 0.5, Lognormal ⇒ 2.0")
for j in 1:8
 plt.subplot(4,2,j)
 #plt.xlabel("time")
 arrival = ""
 service = ""
 μ = 0.0
 if j ∈ (1,3,5,7)
   arrival = "ER0"
   if j ∈ (1,2,3,4)
    μ = 1.0
    if j ∈ (1,2)
     service = "ER1"
    elseif j ∈ (3,4)
     service = "LN1"
    end
   elseif j ∈ (5,6,7,8)
    μ = 2.0
    if j ∈ (5,6)
     service = "ER2"
    elseif j ∈ (7,8)
     service = "LN2"
    end
   end
 elseif j ∈ (2,4,6,8)
   arrival = "LN0"
   if j ∈ (1,2,3,4)
    μ = 1.0
    if j ∈ (1,2)
     service = "ER1"
    elseif j ∈ (3,4)
     service = "LN1"
    end
   elseif j ∈ (5,6,7,8)
    μ = 2.0
    if j ∈ (5,6)
     service = "ER2"
    elseif j ∈ (7,8)
     service = "LN2"
    end
   end
 end
 plt.title("Arrival: dist=$(full_name(arrival)), Service: μ=$μ, dist=$(full_name(service))")

 plt.xlim(0,T)
 #plt.ylim(0,ylim_factor*s)
 queue_length_file_path = "../sum_logs/$queue/$(queue)_$(μ)_$(arrival)_$(service)_$(T)_$(N)_sum_queue_length.txt"
 sojourn_time_file_path = "../sum_logs/$queue/$(queue)_$(μ)_$(arrival)_$(service)_$(T)_$(N)_sum_sojourn_time.txt"
 queue_length_file = open(queue_length_file_path)
 sojourn_time_file = open(sojourn_time_file_path)
 queue_length_file_array = readdlm(queue_length_file)
 sojourn_time_file_array = readdlm(sojourn_time_file)
 close(queue_length_file)
 close(sojourn_time_file)

 X = linspace(0.0,T,T+1)
 #plot_mean_number(queue_length_file_array, plt)
 #plot_ci(queue_length_file_array, plt)
 plot_virtual_sojourn_time(sojourn_time_file_array, plt)
# plot_ci(sojourn_time_file_array, plt)
 plot_sojourn_time_approximation(X, arrival, service, μ, plt)
 plot_sojourn_time_approximation_ρ(X, arrival, service, μ, plt)
 plot_sojourn_time_approximation_GG1(X, arrival, service, μ, plt)
 #plot_arrival_rate_function(X, coeff, plt)
 plt.legend(loc = "lower right", fontsize = 9)
 #plt.tight_layout()
 subplots_adjust(top=0.92,
                 bottom=0.065,
                 left=0.125,
                 right=0.9,
                 hspace=0.305,
                 wspace=0.15)
 plt.savefig("../plots/$(queue)_results.pdf")
end
