type Customer # real customer
  arrival_index::Int64 # this means 'n' of the 'n^th customer'
  remaining_workload::Float64
  arrival_time::Float64 # the time of arrival
  service_beginning_time::Float64 # the time when the service begins
  completion_time::Float64 # the time when the service is completed or when this customer is leaving
  sojourn_time::Float64 # the duration from arrival to leaving
  waiting_time::Float64 # the duration from arrival to the time service begins
end

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
    Λ = (x,s) -> α*(s-x)-(β/γ)*(cos(γ*s)-cos(γ*x))
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
    elseif _control == "DM"
      μ = t -> _TVAS.λ(t) + _param
      M = (x,y) -> _TVAS.Λ(x,y)+(y-x)*_param
      param_letter = "δ"
    end
    string_of_distribution = _string_of_distribution
    new(control, param, param_letter, μ, M, string_of_distribution)
  end
end

type Record
  T::Array{Float64} # time array
  A::Array{Int64}   # number of arrivals
  D::Array{Int64}   # number of departures
  Q::Array{Int64}   # number of customers
  W::Array{Float64} # virtual waiting times (not virtual sojourn time)
  S::Array{Float64} # virtual sojourn times
  file_sim_record::IOStream
  function Record()
    T = Float64[]
    A = Int64[]
    D = Int64[]
    Q = Int64[]
    W = Float64[]
    S = Float64[]
    new(T, A, D, Q, W, S)
  end
end

type TVGG1_queue
  TVAS::Time_Varying_Arrival_Setting
  TVSS::Time_Varying_Service_Setting
  WIP::Array{Customer}
  number_of_customers::Int64
  regular_recording_interval::Float64
  sim_time::Float64
  next_arrival_time::Float64
  next_completion_time::Float64
  next_regular_recording::Float64
  time_index::Int64
  customer_arrival_counter::Int64
  function TVGG1_queue(_TVAS::Time_Varying_Arrival_Setting, _TVSS::Time_Varying_Service_Setting)
    TVAS = _TVAS
    TVSS = _TVSS
    WIP = Customer[]
    number_of_customers = 0
    regular_recording_interval = 0.0
    sim_time = 0.0
    next_arrival_time = typemax(Float64)
    next_completion_time = typemax(Float64)
    next_regular_recording = 0.0
    time_index = 1
    customer_arrival_counter = 0
    new(TVAS, TVSS, WIP, number_of_customers, regular_recording_interval,
        sim_time, next_arrival_time, next_completion_time,
        next_regular_recording, time_index, customer_arrival_counter)
  end
end
