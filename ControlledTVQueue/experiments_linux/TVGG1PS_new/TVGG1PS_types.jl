# Types
type Customer # real customer
  arrival_index::Int64 # this means 'n' of the 'n^th customer'
  remaining_workload::Float64
  arrival_time::Float64 # the time of arrival
  service_beginning_time::Float64 # the time when the service begins
  completion_time::Float64 # the time when the service is completed or when this customer is leaving
  sojourn_time::Float64 # the duration from arrival to leaving
  waiting_time::Float64 # the duration from arrival to the time service begins
end

type Virtual_Customer
  time_index::Int64 # 'n' of the 'n^th virtual customer' == index of time axis. This is used when we record the sojourn time of a virtual customer.
  remaining_workload::Float64
  arrival_time::Float64
  service_beginning_time::Float64
  completion_time::Float64
  sojourn_time::Float64
  waiting_time::Float64
end

type Path
  Status::Array{Customer}
end

type Time_Varying_Arrival_Setting
  α::Float64
  β::Float64
  γ::Float64
  λ::Function # λ(t): arrival rate function
  Λ::Function # Λ(a,b) = ∫λ(s)ds on (a,b]: integrated arrival rate function
  Λ_interval::Function # Λ(a,b) = ∫λ(s)ds on (a,b]: integrated arrival rate function
  string_of_distribution::String # String of base distribution
  base_distribution::Distribution # Nonstationary Non-Poisson process requires its base distribution to generate arrival times
  table::Array{Float64}
  function Time_Varying_Arrival_Setting(_coefficients::Tuple{Float64,Float64,Float64}, _string_of_distribution::String)
    α = _coefficients[1]
    β = _coefficients[2]
    γ = _coefficients[3]
    λ = x -> α + β*sin(γ*x)
    Λ = t -> α*t-(β/γ)*(cos(γ*t)-1)
    Λ_interval = (x,y) -> α*(y-x)-(β/γ)*(cos(γ*y)-cos(γ*x))
    string_of_distribution = _string_of_distribution
    new(α, β, γ, λ, Λ, Λ_interval, string_of_distribution)
  end
end

type Time_Varying_Service_Setting
  control::String
  target::Float64
  string_of_distribution::String # String of workload distribution
  workload_distribution::Distribution # each customer brings its own workload distributed by a certain distribution.
  μ::Function # μ(t): service_rate_function
  M::Function # ∫μ(s)ds on (a,b]: integrated_service_rate_function
  M_interval::Function # ∫μ(s)ds on (a,b]: integrated_service_rate_function
  table::Array{Float64}
  function Time_Varying_Service_Setting(_TVAS::Time_Varying_Arrival_Setting, _control::String, _target::Float64, _string_of_distribution::String)
    control = _control
    target = _target
    string_of_distribution = _string_of_distribution
    new(control, target, string_of_distribution)
  end
end

type Record
  T::Array{Float64} # time array
  A::Array{Int64}   # number of arrivals
  D::Array{Int64}   # number of departures
  Q::Array{Int64}   # number of customers
  W::Array{Float64} # virtual waiting times (not virtual sojourn time)
  S::Array{Float64} # virtual sojourn times
  P::Dict{Float64,Array{Float64}} # Customer path array
  file_sim_record::IOStream
  function Record()
    T = Float64[]
    A = Int64[]
    D = Int64[]
    Q = Int64[]
    W = Float64[]
    S = Float64[]
    P = Dict{Float64, Array{Float64}}()
    new(T, A, D, Q, W, S, P)
  end
end

type TVGG1PS_queue
  TVAS::Time_Varying_Arrival_Setting
  TVSS::Time_Varying_Service_Setting
  WIP::Array{Customer}
  Virtual_WIP::Array{Virtual_Customer}
  number_of_customers::Int64
  number_of_virtual_customers::Int64
  regular_recording_interval::Float64
  sim_time::Float64
  next_arrival_time::Float64
  next_completion_time::Float64
  next_virtual_completion_time::Float64
  next_regular_recording::Float64
  next_completion_index::Int64
  next_virtual_completion_index::Int64
  time_index::Int64
  customer_arrival_counter::Int64
  event::Symbol
  function TVGG1PS_queue(_TVAS::Time_Varying_Arrival_Setting, _TVSS::Time_Varying_Service_Setting)
    TVAS = _TVAS
    TVSS = _TVSS
    WIP = Customer[]
    Virtual_WIP = Virtual_Customer[]
    number_of_customers = 0
    number_of_virtual_customers = 0
    regular_recording_interval = 0.0
    sim_time = 0.0
    next_arrival_time = typemax(Float64)
    next_completion_time = typemax(Float64)
    next_virtual_completion_time = typemax(Float64)
    next_regular_recording = 0.0
    next_completion_index = 0
    next_virtual_completion_index = 0
    time_index = 1
    customer_arrival_counter = 0
    event = :A
    new(TVAS, TVSS, WIP, Virtual_WIP, number_of_customers, number_of_virtual_customers, regular_recording_interval,
        sim_time, next_arrival_time, next_completion_time, next_virtual_completion_time,
        next_regular_recording, next_completion_index, next_virtual_completion_index, time_index, customer_arrival_counter, event)
  end
end
