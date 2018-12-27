include("./functions_distribution.jl")

# Types
mutable struct Customer # real customer
  arrival_index::Int64 # this means 'n' of the 'n^th customer'
  remaining_workload::Float64
  arrival_time::Float64 # the time of arrival
  service_beginning_time::Float64 # the time when the service begins
  completion_time::Float64 # the time when the service is completed or when this customer is leaving
  sojourn_time::Float64 # the duration from arrival to leaving
  waiting_time::Float64 # the duration from arrival to the time service begins
end

mutable struct Time_Varying_Arrival_Setting
  α::Float64
  β::Float64
  γ::Float64
  λ::Function # λ(t): arrival rate function
  Λ::Function # Λ(t) = ∫λ(s)ds on (0,t]: integrated arrival rate function
  Λ_interval::Function # Λ(a,b) = ∫λ(s)ds on (a,b]: integrated arrival rate function
  string_of_distribution::String # String of base distribution
  base_distribution::Distribution # Nonstationary Non-Poisson process requires its base distribution to generate arrival times
  table::Any
  # constructor with sin function
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
  # constructor with general function
  function Time_Varying_Arrival_Setting(_function::Function, _string_of_distribution::String)
    α = 0.0
    β = 0.0
    γ = 0.0
    λ = _function
    Λ = t -> QuadGK.quadgk(λ, 0.0, t)[1]
    Λ_interval = (x,y) -> QuadGK.quadgk(λ, x, y)[1]
    string_of_distribution = _string_of_distribution
    new(α, β, γ, λ, Λ, Λ_interval, string_of_distribution)
  end
end

mutable struct Time_Varying_Service_Setting
  control::String
  target::Float64
  string_of_distribution::String # String of workload distribution
  workload_distribution::Distribution # each customer brings its own workload distributed by a certain distribution.
  μ::Function # μ(t): service_rate_function
  M::Function # ∫μ(s)ds on (a,b]: integrated_service_rate_function
  M_interval::Function # ∫μ(s)ds on (a,b]: integrated_service_rate_function
  table::Any
  function Time_Varying_Service_Setting(_TVAS::Time_Varying_Arrival_Setting, _control::String, _target::Float64, _string_of_distribution::String)
    control = _control
    target = _target
    string_of_distribution = _string_of_distribution
    new(control, target, string_of_distribution)
  end
end

mutable struct TVGG1PS_queue
  TVAS::Time_Varying_Arrival_Setting
  TVSS::Time_Varying_Service_Setting
  WIP::Array{Customer}
  number_of_customers::Int64
  regular_recording_interval::Float64
  sim_time::Float64
  next_arrival_time::Float64
  next_completion_time::Float64
  next_regular_recording::Float64
  next_completion_index::Int64
  time_index::Int64
  customer_arrival_counter::Int64
  event::Symbol
  size_customer_pool::Int64
  function TVGG1PS_queue(_TVAS::Time_Varying_Arrival_Setting, _TVSS::Time_Varying_Service_Setting)
    TVAS = _TVAS
    TVSS = _TVSS
    WIP = Customer[]
    number_of_customers = 0
    regular_recording_interval = 0.0
    sim_time = 0.0
    next_arrival_time = typemax(Float64)
    next_completion_time = typemax(Float64)
    next_regular_recording = 0.0
    next_completion_index = 0
    time_index = 1
    customer_arrival_counter = 0
    event = :A
    size_customer_pool = 0
    new(TVAS, TVSS, WIP, number_of_customers, regular_recording_interval,
        sim_time, next_arrival_time, next_completion_time,
        next_regular_recording, next_completion_index, time_index, customer_arrival_counter, event, size_customer_pool)
  end
end

mutable struct TVGG1PS_queue_path
  WIP::Array{Customer}
  number_of_customers::Int64
  regular_recording_interval::Float64
  sim_time::Float64
  next_arrival_time::Float64
  next_completion_time::Float64
  next_regular_recording::Float64
  next_completion_index::Int64
  time_index::Int64
  customer_arrival_counter::Int64
  event::Symbol
  size_customer_pool::Int64
  function TVGG1PS_queue_path(system::TVGG1PS_queue)
    WIP = deepcopy(system.WIP)
    number_of_customers = system.number_of_customers
    regular_recording_interval = system.regular_recording_interval
    sim_time = system.sim_time
    next_arrival_time = system.next_arrival_time
    next_completion_time = system.next_completion_time
    next_regular_recording = system.next_regular_recording
    next_completion_index = system.next_completion_index
    time_index = system.time_index
    customer_arrival_counter = system.customer_arrival_counter
    event = system.event
    size_customer_pool = system.size_customer_pool
    new(WIP, number_of_customers, regular_recording_interval,
        sim_time, next_arrival_time, next_completion_time,
        next_regular_recording, next_completion_index, time_index, customer_arrival_counter, event, size_customer_pool)
  end
end

mutable struct Record
  T::Array{Float64} # time array
  A::Array{Int64}   # number of arrivals
  D::Array{Int64}   # number of departures
  Q::Array{Int64}   # number of customers
  W::Array{Float64} # virtual waiting times (not virtual sojourn time)
  S::Array{Float64} # virtual sojourn times
  P::Dict{Float64,TVGG1PS_queue_path} # path array
  file_sim_record::IOStream
  function Record()
    T = Float64[]
    A = Int64[]
    D = Int64[]
    Q = Int64[]
    W = Float64[]
    S = Float64[]
    P = Dict{Float64, TVGG1PS_queue_path}()
    new(T, A, D, Q, W, S, P)
  end
end
