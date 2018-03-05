function table_for_any_periodic_J(f::Function, T::Float64)  # standard_J: λ(t)=1+βsin(t)
    temp = Float64[]
    X = 0.0:T/10000:T
    for x in X push!(temp,f(x)) end
    f_u = maximum(temp)
    f_l = minimum(temp)
    average = QuadGK.quadgk(f,0.0,T)[1]/T
    F(t) = QuadGK.quadgk(f,0.0,t)[1]
    ϵ = max(1,1/f_u)/100
    ρ = f_u/f_l
    η = ϵ/(1+ρ)
    δ = (f_u*ϵ)/(1+ρ)
    nx = (T*(1+ρ))/ϵ
    ny = F(T)/δ
    y = 0.0:δ:F(T)

    X = 0.0:η:T
    a = Float64[]
    for x in X
        push!(a,F(x))
    end

    b = Float64[]
    i,j = 1,1
    while j < nx + 1 && i < ny + 1
        if y[i] > a[j]
            j += 1
        else
            push!(b, j)
            i += 1
        end
    end
    push!(b,average)    # b[end-3]
    push!(b,T)        # b[end-2]
    push!(b,δ)        # b[end-1]
    push!(b,η)        # b[end]
    return b
end

function table_for_sinusoidal_J(coeff::Tuple{Float64,Float64,Float64})  # standard_J: λ(t)=1+βsin(t)
    α, β, γ = coeff
    λ_bar = α
    T = (2*π)/γ
    λ(t) = α+β*sin(γ*t)
    Λ(t) = α*t-(β/γ)*(cos(γ*t)-1)
    λ_u, λ_l = α+β, α-β
    ϵ = max(1,1/λ_u)/100
    ρ = λ_u/λ_l
    η = ϵ/(1+ρ)
    δ = (λ_u*ϵ)/(1+ρ)
    nx = (T*(1+ρ))/ϵ
    ny = Λ(T)/δ
    y = 0.0:δ:Λ(T)
    a = Λ(0.0:η:T)
    b = Float64[]
    i,j = 1,1
    while j < nx + 1 && i < ny + 1
        if y[i] > a[j]
            j += 1
        else
            push!(b, j)
            i += 1
        end
    end
    push!(b,λ_bar)    # b[end-3]
    push!(b,T)        # b[end-2]
    push!(b,δ)        # b[end-1]
    push!(b,η)        # b[end]
    return b
end

function J(t::Float64, table::Array{Float64})
    index = Int64(floor(mod(t,(table[end-3]*table[end-2]))/table[end-1]))
    if index != 0
        return floor(t/(table[end-3]*table[end-2]))*table[end-2] + table[index]*table[end]
    else
        return floor(t/(table[end-3]*table[end-2]))*table[end-2]
    end
end

function inverse_integrated_rate_function(Λ::Function, s::Float64, val::Float64, table::Array{Float64})
    return J(val+Λ(s),table)
end

function inverse_cdf(f::Function, val::Float64) # find inf{t>0: F(t) > val}, f: pdf
  s = 0.0 # integration start point
  t = 0.0 # integration end point
  temp = 0.0 # temporal integration value (s,t)
  sum = 0.0 # whole integration value (0.0,t)
  while sum < val
    temp = QuadGK.quadgk(f,s,t)[1] # compute ∫ s to t
    s = t # change the value of s to t
    t += 0.001 # increase the value of t by 0.001
    sum += temp # add temp to sum
  end
  return t
end

function string_to_dist(Distribution::String) # All distributions have mean 1.0.
  if Distribution == "Pareto"
    return Pareto(1+sqrt(2),2-sqrt(2))
  elseif Distribution == "Lognormal"
    return LogNormal(-log(2)/2,sqrt(log(2)))
  elseif Distribution == "EXP"
    return Exponential(1.0)
  elseif Distribution == "E2"
    return Erlang(2,1/2)
  elseif Distribution == "H2"
    return Exponential(1.0) # this is temporal
  end
end

function set_distribution(TVS::Any)
  if typeof(TVS) == Time_Varying_Arrival_Setting
    TVS.base_distribution = string_to_dist(TVS.string_of_distribution)
  elseif typeof(TVS) == Time_Varying_Service_Setting
    TVS.workload_distribution = string_to_dist(TVS.string_of_distribution)
  end
end

function set_service_rate_function(TVAS::Any, TVSS::Any)
  scv_arrival = 0.0
  scv_service = 0.0
  if TVAS.string_of_distribution == "EXP"
    scv_arrival = 1.0
  elseif TVAS.string_of_distribution == "H2"
    scv_arrival = 4.0
  elseif TVAS.string_of_distribution == "E2"
    scv_arrival = 0.5
  end

  if TVSS.string_of_distribution == "EXP"
    scv_service = 1.0
  elseif TVSS.string_of_distribution == "H2"
    scv_service = 4.0
  elseif TVSS.string_of_distribution == "E2"
    scv_service = 0.5
  end

  s = TVSS.target
  V = (scv_arrival+scv_service)/2
  λ = TVAS.λ
  if TVSS.control == "SR"
    TVSS.μ = t -> λ(t)/2-(1/(2*s))+sqrt(((s*λ(t)+1)^2)-4*s*(1-V)*λ(t))/(2*s)
    TVSS.M = t -> QuadGK.quadgk(TVSS.μ, 0.0, t)[1]
    TVSS.M_interval = (x,y) -> QuadGK.quadgk(TVSS.μ, x, y)[1]
  elseif TVSS.control == "PD"
    TVSS.μ = t -> λ(t) + (V/s)
    TVSS.M = t -> TVAS.Λ(t)+t*(V/s)
    TVSS.M_interval = (x,y) -> TVAS.Λ_interval(x,y)+(y-x)*(V/s)
  end
end

function set_tables(TVAS::Any, TVSS::Any)
    TVAS.table = table_for_sinusoidal_J((TVAS.α,TVAS.β,TVAS.γ))
    TVSS.table = table_for_any_periodic_J(TVSS.μ, 2*π/TVAS.γ)
end

function generate_Hyperexponential(p1::Float64, p2::Float64, θ1::Float64, θ2::Float64)
  return rand() < p1 ? rand(Exponential(θ1)) : rand(Exponential(θ2)) # with prob. p1, return Exponential(1/λ1), with prob. p2, returen Exp(1/λ2)
end

function generate_NHNP(TVAS::Time_Varying_Arrival_Setting, T::Float64)
  if TVAS.string_of_distribution != "H2"
    g_e = t -> 1-cdf(TVAS.base_distribution,t)/mean(TVAS.base_distribution)
    S = inverse_cdf(g_e , rand())
  	V = [ inverse_integrated_rate_function(TVAS.Λ, 0.0, S, TVAS.table) ]
    n = 1
    while V[n] < T
      push!(V, inverse_integrated_rate_function(TVAS.Λ, V[n], rand(TVAS.base_distribution), TVAS.table) )
      n += 1
    end
    return V
  elseif TVAS.string_of_distribution == "H2"
    # Hyperexponential parameters
    p1 = (5+sqrt(15))/10
    p2 = 1-p1
    θ1 = 1/(2*p1)
    θ2 = 1/(2*p2)
    g = t -> p1*(1/θ1)*exp(-(1/θ1)*t)+p2*(1/θ2)*exp(-(1/θ2)*t) # Hyperexponential pdf
    G = t -> QuadGK.quadgk(g, 0.0, t)[1] # Hyperexponential cdf
    τ = 1.0 # Hyeperexponential mean
    g_e = t -> (1-G(t))/τ # Hyperexponential equilibrium pdf
    S = inverse_cdf(g_e, rand())
    V = [ inverse_integrated_rate_function(TVAS.Λ, 0.0, S, TVAS.table) ]
    n = 1
    while V[n] < T
      S = generate_Hyperexponential(p1,p2,θ1,θ2)
      push!(V, inverse_integrated_rate_function(TVAS.Λ, V[n], S, TVAS.table) )
      n += 1
    end
    return V
  end
end

function generate_customer_pool(TVAS::Time_Varying_Arrival_Setting, TVSS::Time_Varying_Service_Setting, T::Float64)
  if TVSS.string_of_distribution != "H2"
    arrival_times = generate_NHNP(TVAS, T)
    Customers = Customer[]
    for i in 1:length(arrival_times)
      push!(Customers, Customer(i, rand(TVSS.workload_distribution), arrival_times[i], 0.0, typemax(Float64), 0.0, 0.0))
    end
    return Customers
  else TVSS.string_of_distribution == "H2"
    p1 = (5+sqrt(15))/10
    p2 = 1-p1
    θ1 = 1/(2*p1)
    θ2 = 1/(2*p2)
    arrival_times = generate_NHNP(TVAS, T)
    Customers = Customer[]
    for i in 1:length(arrival_times)
      push!(Customers, Customer(i, generate_Hyperexponential(p1,p2,θ1,θ2), arrival_times[i], 0.0, typemax(Float64), 0.0, 0.0))
    end
    return Customers
  end
end

function next_event(system::TVGG1_queue, customer_pool::Array{Customer}, record::Record)
  if system.next_arrival_time == min(system.next_arrival_time, system.next_completion_time, system.next_regular_recording)
    #println("Event: NEW ARRIVAL")
    current_time = system.next_arrival_time
    # reduce workloads for each customer & virtual customer in system
    service_amount = 0.0
    if system.number_of_customers > 0
      service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      system.WIP[1].remaining_workload -= service_amount
    end

    # insert a new customer in the system
    system.customer_arrival_counter += 1
    system.number_of_customers += 1
    push!(system.WIP, customer_pool[system.customer_arrival_counter])

    # update simulational time
    system.sim_time = current_time

    # update next arrival time
    system.next_arrival_time = customer_pool[system.customer_arrival_counter+1].arrival_time

    # update next completion information
    if system.number_of_customers == 1
      system.WIP[1].service_beginning_time = system.WIP[1].arrival_time
      system.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, system.WIP[1].remaining_workload, system.TVSS.table)
    end
  elseif system.next_completion_time == min(system.next_arrival_time, system.next_completion_time, system.next_regular_recording)
    current_time = system.next_completion_time

    # save completed Customer time information
    system.WIP[1].completion_time = current_time
    system.WIP[1].sojourn_time = current_time - system.WIP[1].arrival_time

    # remove the completing customer from the system
    shift!(system.WIP)
    system.number_of_customers -= 1

    # update simulational time
    system.sim_time = current_time

    # update next completion information & next customer's service beginning time
    if system.number_of_customers > 0
      system.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, system.WIP[1].remaining_workload, system.TVSS.table)
      system.WIP[1].service_beginning_time = system.sim_time
      system.WIP[1].waiting_time = system.sim_time - system.WIP[1].arrival_time
    else
      system.next_completion_time = typemax(Float64)
    end
  elseif system.next_regular_recording == min(system.next_arrival_time, system.next_completion_time, system.next_regular_recording)
    current_time = system.next_regular_recording

    # reduce workloads for a customer
    if system.number_of_customers > 0
      service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      system.WIP[1].remaining_workload -= service_amount
    end

    # update simulational time & increase time index
    system.sim_time = current_time
    system.time_index += 1

    # update next regular recording time
    system.next_regular_recording += system.regular_recording_interval

    # record the number of customer, arrival process, temporal sojourn time
    push!(record.Q, system.number_of_customers)
    push!(record.A, system.customer_arrival_counter)
  end
end

function run_to_end(system::TVGG1_queue, customer_pool::Array{Customer}, record::Record, T::Float64, warm_up_time::Float64)
  while system.sim_time < T
    next_event(system, customer_pool, record)
  end
end

function record_virtual_waiting_times(customer_pool::Array{Customer}, record::Record)
  for i in 1:length(record.T)
    t = record.T[i] # time t
    if record.A[i] == 0
      push!(record.W, 0.0)
    else
      k = record.A[i] # customer index at time t
      push!(record.W, max(0.0, customer_pool[k].sojourn_time - (t - customer_pool[k].arrival_time)) )
    end
  end
end

function record_virtual_sojourn_times(TVSS::Time_Varying_Service_Setting, record::Record)
  if TVSS.string_of_distribution != "H2"
    for i in 1:length(record.T)
      push!(record.S , record.W[i] + inverse_integrated_rate_function(TVSS.M, record.T[i] + record.W[i], rand(TVSS.workload_distribution), TVSS.table) - (record.T[i] + record.W[i]) )
    end
  elseif TVSS.string_of_distribution == "H2"
    p1 = (5+sqrt(15))/10
    p2 = 1-p1
    θ1 = 1/(2*p1)
    θ2 = 1/(2*p2)
    for i in 1:length(record.T)
      push!(record.S , record.W[i] + inverse_integrated_rate_function(TVSS.M, record.T[i] + record.W[i], generate_Hyperexponential(p1,p2,θ1,θ2), TVSS.table) - (record.T[i] + record.W[i]) )
    end
  end
end

function do_experiment(queue::String, control::String, target::Float64, arrival::String,
                        service::String, coeff::Tuple, T::Float64, N::Int64, record::Record)
  TVAS = Time_Varying_Arrival_Setting(coeff, arrival)
  TVSS = Time_Varying_Service_Setting(TVAS, control, target, service)
  set_distribution(TVAS)
  set_distribution(TVSS)
  set_service_rate_function(TVAS, TVSS)
  set_tables(TVAS,TVSS)
  file_num_in_queue = open("../../../logs/queue_length_$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)).txt" , "w")
  file_virtual_sojourn_time = open("../../../logs/sojourn_time_$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)).txt" , "w")
  regular_recording_interval = T/1000
  t = 0.0
  while t <= T
		push!(record.T, t)
    t += regular_recording_interval
  end
  writedlm(file_num_in_queue, transpose(record.T)) # write time-axis
  writedlm(file_virtual_sojourn_time, transpose(record.T)) # write time-axis

  for n in 1:N
    println("Replication $n")
	record.Q = Int64[]
    record.W = Float64[]
    record.S = Float64[]
    record.A = Int64[]
    system = TVGG1_queue(TVAS, TVSS)
    customer_pool = generate_customer_pool(TVAS, TVSS, T*1.5)
    system.regular_recording_interval = T/1000
    system.next_arrival_time = customer_pool[1].arrival_time
    run_to_end(system, customer_pool, record, T*1.2, 0.0)
    record_virtual_waiting_times(customer_pool, record)
    record_virtual_sojourn_times(TVSS, record)
    writedlm(file_num_in_queue, transpose(record.Q)) # write record Q(t)
    writedlm(file_virtual_sojourn_time, transpose(record.S)) # write record S(t)
	end
  close(file_num_in_queue)
  close(file_virtual_sojourn_time)
end
