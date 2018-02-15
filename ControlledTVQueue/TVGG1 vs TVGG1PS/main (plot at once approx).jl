using PyPlot, Distributions, QuadGK

type Time_Varying_Arrival_Setting
  α::Float64
  β::Float64
  γ::Float64
  λ::Function # λ(t): arrival rate function
  Λ::Function # Λ(a,b) = ∫λ(s)ds on (a,b]: integrated arrival rate function
  string_of_distribution::String # String of base distribution
  base_distribution::Distribution # Nonstationary Non-Poisson process requires its base distribution to generate arrival times
  function Time_Varying_Arrival_Setting(_coefficients::Tuple{Float64,Float64,Float64}, _string_of_distribution::String)
    α = _coefficients[1]
    β = _coefficients[2]
    γ = _coefficients[3]
    λ = x -> α + β*sin(γ*x)
    Λ = (x,s) -> α*(s-x)-(β/γ)*(cos(γ*x)-cos(γ*s))
    string_of_distribution = _string_of_distribution
    new(α, β, γ, λ, Λ, string_of_distribution)
  end
end

type Time_Varying_Service_Setting
  control::String
  param::Float64
  param_letter::String
  μ::Function # μ(t): service_rate_function
  M::Function # ∫μ(s)ds on (a,b]: integrated_service_rate_function
  string_of_distribution::String # String of workload distribution
  workload_distribution::Distribution # each customer brings its own workload distributed by a certain distribution.
  function Time_Varying_Service_Setting(_TVAS::Time_Varying_Arrival_Setting, _control::String, _param::Float64, _string_of_distribution::String)
    control = _control
    param = _param
    if _control == "RM"
      μ = t -> _TVAS.λ(t)/_param
      M = (x,y) -> (1/_param)*_TVAS.Λ(x,y)
      param_letter = "ρ"
    elseif _control == "SR"
      μ = t -> _TVAS.λ(t) + _param*sqrt(_TVAS.λ(t))
      M = (x,y) -> QuadGK.quadgk(μ, x, y)[1]
      param_letter = "ξ"
    elseif _control == "PSA"
      μ = t -> _TVAS.λ(t) + (_TVAS.λ(t)/2)*(sqrt((_TVAS.λ(t)+_param)/_TVAS.λ(t))-1)
      M = (x,y) -> QuadGK.quadgk(μ, x, y)[1]
      param_letter = "ζ"
    end
    string_of_distribution = _string_of_distribution
    new(control, param, param_letter, μ, M, string_of_distribution)
  end
end

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

function plot_arrival_rate_function(x::Any, TVAS::Time_Varying_Arrival_Setting, plt::Module)
  λ = TVAS.λ(x)
  plt.plot(x,λ,color="black",linestyle="--", label="arrival rate")
end

function plot_apprx_mean_number(x::Any, TVAS::Time_Varying_Arrival_Setting, TVSS::Time_Varying_Service_Setting, plt::Module, queue::String)
  Q = Float64[]
  if queue == "TVGG1"
    ρ = t -> TVAS.λ(t)/TVSS.μ(t)
    for t in x
      push!(Q , ρ(t)+((scv_of_distribution(TVAS.string_of_distribution)+scv_of_distribution(TVSS.string_of_distribution))/2)*((ρ(t)^2)/(1-ρ(t))) )
    end
  elseif queue == "TVGG1PS"
    ρ = t -> TVAS.λ(t)/TVSS.μ(t)
    for t in x
      push!(Q , ((scv_of_distribution(TVAS.string_of_distribution)+scv_of_distribution(TVSS.string_of_distribution))/(1+scv_of_distribution(TVSS.string_of_distribution)))*(ρ(t)/(1-ρ(t))) )
    end
  end
  plt.plot(x,Q,color="green", label="E[Q(t)]")
end

function plot_apprx_mean_sojourn_time(x::Any, TVAS::Time_Varying_Arrival_Setting, TVSS::Time_Varying_Service_Setting, plt::Module, queue::String)
  S = Float64[]
  if queue == "TVGG1"
    ρ = t -> TVAS.λ(t)/TVSS.μ(t)
    for t in x
      push!(S , (1/TVAS.λ(t))*(ρ(t)+((scv_of_distribution(TVAS.string_of_distribution)+scv_of_distribution(TVSS.string_of_distribution))/2)*((ρ(t)^2)/(1-ρ(t)))) )
    end
  elseif queue == "TVGG1PS"
    ρ = t -> TVAS.λ(t)/TVSS.μ(t)
    #β = (1/TVSS.μ(t)) # mean job size
    for t in x
      push!(S , ((scv_of_distribution(TVAS.string_of_distribution)+scv_of_distribution(TVSS.string_of_distribution))/(1+scv_of_distribution(TVSS.string_of_distribution)))*((1/TVSS.μ(t))/(1-ρ(t))) )
    end
  end
  plt.plot(x,S,color="red", label="E[S(t)]")
end

function plot_apprx_mean_number(x::Any, TVAS::Time_Varying_Arrival_Setting, TVSS::Time_Varying_Service_Setting, plt::Module, _label::String, _color::String, queue::String)
  Q = Float64[]
  if queue == "TVGG1"
    ρ = t -> TVAS.λ(t)/TVSS.μ(t)
    for t in x
      push!(Q , ρ(t)+((scv_of_distribution(TVAS.string_of_distribution)+scv_of_distribution(TVSS.string_of_distribution))/2)*((ρ(t)^2)/(1-ρ(t))) )
    end
  elseif queue == "TVGG1PS"
    ρ = t -> TVAS.λ(t)/TVSS.μ(t)
    for t in x
      push!(Q , ((scv_of_distribution(TVAS.string_of_distribution)+scv_of_distribution(TVSS.string_of_distribution))/(1+scv_of_distribution(TVSS.string_of_distribution)))*(ρ(t)/(1-ρ(t))) )
    end
  end
  plt.plot(x,Q,color=_color, label=_label)
end

function plot_apprx_mean_sojourn_time(x::Any, TVAS::Time_Varying_Arrival_Setting, TVSS::Time_Varying_Service_Setting, plt::Module, _label::String, _color::String, queue::String)
  S = Float64[]
  if queue == "TVGG1"
    ρ = t -> TVAS.λ(t)/TVSS.μ(t)
    for t in x
      push!(S , (1/TVAS.λ(t))*(ρ(t)+((scv_of_distribution(TVAS.string_of_distribution)+scv_of_distribution(TVSS.string_of_distribution))/2)*((ρ(t)^2)/(1-ρ(t)))) )
    end
  elseif queue == "TVGG1PS"
    ρ = t -> TVAS.λ(t)/TVSS.μ(t)
    #β = 1.0 # mean job size
    for t in x
      push!(S , ((scv_of_distribution(TVAS.string_of_distribution)+scv_of_distribution(TVSS.string_of_distribution))/(1+scv_of_distribution(TVSS.string_of_distribution)))*((1/TVSS.μ(t))/(1-ρ(t))) )
    end
  end
  plt.plot(x,S,color=_color, label=_label)
end

function scv_of_distribution(string_of_distribution::String)
  if string_of_distribution == "Exponential"
    return 1.0
  elseif string_of_distribution == "Hyperexponential"
    return 4.0
  elseif string_of_distribution == "Erlang"
    return 0.5
  end
end

# ALL AT ONCE
queue_set = ("TVGG1","TVGG1PS")
queue_type_set = (L"$G_t/G_t/1$", L"$G_t/G_t/1/PS$")

control_set = ("RM", "SR", "PSA")
param_letter_set = ("ρ", "ξ", "ζ")
param_set = (0.8,0.2,1.0)
distribution_set = [("Exponential","Exponential"), ("Hyperexponential","Hyperexponential"), ("Hyperexponential","Erlang"), ("Erlang","Erlang")]
short_distribution_set = [("EXP","EXP"), ("H2","H2"), ("H2","E2"), ("E2","E2")]
gamma_set = (0.001, 0.01, 0.1, 1.0 , 10.0)
time_set = (20000.0, 2000.0, 2000.0, 2000.0, 2000.0)
colors = ("green","red","blue","magenta")

ylim_set_1 = Any[]
push!(ylim_set_1, [(5.0,5.0,5.0,5.0,5.0),(18.0,18.0,18.0,18.0,18.0),(12.0,12.0,12.0,12.0,12.0),(4.0,4.0,4.0,4.0,4.0)])
push!(ylim_set_1, [(7.0,7.0,7.0,7.0,7.0),(20.0,20.0,20.0,20.0,20.0),(15.0,15.0,15.0,15.0,15.0),(4.0,4.0,4.0,4.0,4.0)])
push!(ylim_set_1, [(6.0,6.0,6.0,6.0,6.0),(22.0,22.0,22.0,22.0,22.0),(16.0,16.0,16.0,16.0,16.0),(4.0,4.0,4.0,4.0,4.0)])
ylim_set_2=[(14.0,14.0,14.0,14.0,14.0),(20.0,20.0,20.0,20.0,20.0),(21.0,21.0,21.0,21.0,21.0)]
ylim_set_3=[(16.0,16.0,16.0,16.0,16.0),(20.0,20.0,20.0,20.0,20.0),(20.0,20.0,20.0,20.0,20.0)]

for c in 1:length(control_set)
 param_letter = param_letter_set[c]
 param = param_set[c]
 control = control_set[c]
 for d in 1:length(distribution_set)
   # apprx Plot type 1 (mean number and wait)
   plt = PyPlot
   plt.figure(figsize=(10.5,15))
   arrival = distribution_set[d][1]
   service = distribution_set[d][2]
   short_arrival = short_distribution_set[d][1]
   short_service = short_distribution_set[d][2]

   queue = queue_set[1]
   queue_type = queue_type_set[1]
   for i in 1:length(gamma_set)
     plt.subplot(5,2,(2*i)-1)
     coeff = (1.0, 0.2, gamma_set[i])
     time = time_set[i]
     TVAS = Time_Varying_Arrival_Setting(coeff, arrival)
     TVSS = Time_Varying_Service_Setting(TVAS, control, param, service)
     plt.title("$queue_type, $control ($param_letter = $param), $short_arrival/$short_service, γ = $(coeff[3])")
     plt.xlabel("time")
     plt.xlim(0,time)
     plt.ylim(0,ylim_set_1[c][d][i])
     x = linspace(0.0, time, 100000)
     plot_apprx_mean_number(x, TVAS, TVSS, plt, queue)
     plot_apprx_mean_sojourn_time(x, TVAS, TVSS, plt, queue)
     plot_arrival_rate_function(x, TVAS, plt)
     plt.legend(loc = "right", fontsize = 8)
   end

   queue = queue_set[2]
   queue_type = queue_type_set[2]
   for i in 1:length(gamma_set)
     plt.subplot(5,2,(2*i))
     coeff = (1.0, 0.2, gamma_set[i])
     time = time_set[i]
     TVAS = Time_Varying_Arrival_Setting(coeff, arrival)
     TVSS = Time_Varying_Service_Setting(TVAS, control, param, service)
     plt.title("$queue_type, $control ($param_letter = $param), $short_arrival/$short_service, γ = $(coeff[3])")
     plt.xlabel("time")
     plt.xlim(0,time)
     plt.ylim(0,ylim_set_1[c][d][i])
     x = linspace(0.0, time, 100000)
     plot_apprx_mean_number(x, TVAS, TVSS, plt, queue)
     plot_apprx_mean_sojourn_time(x, TVAS, TVSS, plt, queue)
     plot_arrival_rate_function(x, TVAS, plt)
     plt.legend(loc = "right", fontsize = 8)
   end

   plt.tight_layout()
   plt.savefig("./plot type 1 (apprx mean number and sojourn time) (TVGG1 vs TVGG1PS, $control, $arrival, $service).pdf")
   close()
 end

 # apprx Plot type 3 (mean number for general distributions)
 plt = PyPlot
 plt.figure(figsize=(10.5,15))
 queue = queue_set[1]
 queue_type = queue_type_set[1]
 for i in 1:length(gamma_set)
   γ = gamma_set[i]
   time = time_set[i]
   x = linspace(0.0, time, 10000)
   ## plotting
   plt.subplot(5,2,(2*i)-1)
   coeff = (1.0, 0.2, γ)
   plt.title("$queue_type, $control ($param_letter = $param), γ = $(coeff[3])")
   plt.xlabel("time")
   plt.xlim(0,time)
   plt.ylim(0,ylim_set_2[c][i])
   distribution_set_2 = [("Exponential","Exponential"),("Hyperexponential","Hyperexponential"),("Hyperexponential","Erlang"),("Erlang","Erlang")]
   labels = ("E[Q(t)] (EXP/EXP)","E[Q(t)] (H2/H2)","E[Q(t)] (H2/E2)","E[Q(t)] (E2/E2)")
   for j in 1:length(distribution_set_2)
     TVAS = Time_Varying_Arrival_Setting(coeff, distribution_set_2[j][1])
     TVSS = Time_Varying_Service_Setting(TVAS, control, param, distribution_set_2[j][2])
     plot_apprx_mean_number(x, TVAS, TVSS, plt, labels[j], colors[j], queue)
   end
   TVAS = Time_Varying_Arrival_Setting(coeff, distribution_set_2[1][1])
   plot_arrival_rate_function(x, TVAS, plt)
   plt.legend(loc = "upper right", fontsize = 8)
 end

 queue = queue_set[2]
 queue_type = queue_type_set[2]
 for i in 1:length(gamma_set)
   γ = gamma_set[i]
   time = time_set[i]
   x = linspace(0.0, time, 10000)
   ## plotting
   plt.subplot(5,2,(2*i))
   coeff = (1.0, 0.2, γ)
   plt.title("$queue_type, $control ($param_letter = $param), γ = $(coeff[3])")
   plt.xlabel("time")
   plt.xlim(0,time)
   plt.ylim(0,ylim_set_2[c][i])
   distribution_set_2 = [("Exponential","Exponential"),("Hyperexponential","Hyperexponential"),("Hyperexponential","Erlang"),("Erlang","Erlang")]
   labels = ("E[Q(t)] (EXP/EXP)","E[Q(t)] (H2/H2)","E[Q(t)] (H2/E2)","E[Q(t)] (E2/E2)")
   for j in 1:length(distribution_set_2)
     TVAS = Time_Varying_Arrival_Setting(coeff, distribution_set_2[j][1])
     TVSS = Time_Varying_Service_Setting(TVAS, control, param, distribution_set_2[j][2])
     plot_apprx_mean_number(x, TVAS, TVSS, plt, labels[j], colors[j], queue)
   end
   TVAS = Time_Varying_Arrival_Setting(coeff, distribution_set_2[1][1])
   plot_arrival_rate_function(x, TVAS, plt)
   plt.legend(loc = "upper right", fontsize = 8)
 end

 plt.tight_layout()
 plt.savefig("./plot type 3 (apprx mean number) (TVGG1 vs TVGG1PS, $control, general distributions).pdf")
 close()

 # apprx Plot type 4 (mean sojourn time for general distributions)
 plt.figure(figsize=(10.5,15))
 queue = queue_set[1]
 queue_type = queue_type_set[1]
 for i in 1:length(gamma_set)
   γ = gamma_set[i]
   time = time_set[i]
   x = linspace(0.0, time, 10000)
   ## plotting
   plt.subplot(5,2,(2*i)-1)
   coeff = (1.0, 0.2, γ)
   plt.title("$queue_type, $control ($param_letter = $param), γ = $(coeff[3])")
   plt.xlabel("time")
   plt.xlim(0,time)
   plt.ylim(0,ylim_set_3[c][i])
   distribution_set_2 = [("Exponential","Exponential"),("Hyperexponential","Hyperexponential"),("Hyperexponential","Erlang"),("Erlang","Erlang")]
   labels = ("E[S(t)] (EXP/EXP)","E[S(t)] (H2/H2)","E[S(t)] (H2/E2)","E[S(t)] (E2/E2)")
   for j in 1:length(distribution_set_2)
     TVAS = Time_Varying_Arrival_Setting(coeff, distribution_set_2[j][1])
     TVSS = Time_Varying_Service_Setting(TVAS, control, param, distribution_set_2[j][2])
     plot_apprx_mean_sojourn_time(x, TVAS, TVSS, plt, labels[j], colors[j], queue)
   end
   TVAS = Time_Varying_Arrival_Setting(coeff, distribution_set_2[1][1])
   plot_arrival_rate_function(x, TVAS, plt)
   plt.legend(loc = "upper right", fontsize = 8)
 end

 queue = queue_set[2]
 queue_type = queue_type_set[2]
 for i in 1:length(gamma_set)
   γ = gamma_set[i]
   time = time_set[i]
   x = linspace(0.0, time, 10000)
   ## plotting
   plt.subplot(5,2,(2*i))
   coeff = (1.0, 0.2, γ)
   plt.title("$queue_type, $control ($param_letter = $param), γ = $(coeff[3])")
   plt.xlabel("time")
   plt.xlim(0,time)
   plt.ylim(0,ylim_set_3[c][i])
   distribution_set_2 = [("Exponential","Exponential"),("Hyperexponential","Hyperexponential"),("Hyperexponential","Erlang"),("Erlang","Erlang")]
   labels = ("E[S(t)] (EXP/EXP)","E[S(t)] (H2/H2)","E[S(t)] (H2/E2)","E[S(t)] (E2/E2)")
   for j in 1:length(distribution_set_2)
     TVAS = Time_Varying_Arrival_Setting(coeff, distribution_set_2[j][1])
     TVSS = Time_Varying_Service_Setting(TVAS, control, param, distribution_set_2[j][2])
     plot_apprx_mean_sojourn_time(x, TVAS, TVSS, plt, labels[j], colors[j], queue)
   end
   TVAS = Time_Varying_Arrival_Setting(coeff, distribution_set_2[1][1])
   plot_arrival_rate_function(x, TVAS, plt)
   plt.legend(loc = "upper right", fontsize = 8)
 end

 plt.tight_layout()
 plt.savefig("./plot type 4 (apprx mean sojourn time) (TVGG1 vs TVGG1PS, $control, general distributions).pdf")
 close()
end
