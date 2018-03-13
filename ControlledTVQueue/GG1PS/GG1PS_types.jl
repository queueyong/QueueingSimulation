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

type Arrival_Setting
  string_of_distribution::String # String of interarrival distribution
  λ::Float64 # arrival rate
  interarrival_distribution::Distribution # Distributions.distribution
  function Arrival_Setting(dist::String)
    string_of_distribution = dist
    new(string_of_distribution)
  end
end

type Service_Setting
  string_of_distribution::String # String of workload distribution
  μ::Float64 # service rate
  workload_distribution::Distribution # each customer brings its own workload distributed by a certain distribution.
  function Service_Setting(dist::String, _μ::Float64)
    string_of_distribution = dist
    μ = _μ
    new(string_of_distribution, μ)
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

type GG1PS_queue
  AS::Arrival_Setting
  SS::Service_Setting
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
  function GG1PS_queue(_AS::Arrival_Setting, _SS::Service_Setting)
    AS = _AS
    SS = _SS
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
    new(AS, SS, WIP, Virtual_WIP, number_of_customers, number_of_virtual_customers, regular_recording_interval,
        sim_time, next_arrival_time, next_completion_time, next_virtual_completion_time,
        next_regular_recording, next_completion_index, next_virtual_completion_index, time_index, customer_arrival_counter)
  end
end
