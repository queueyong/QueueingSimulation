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
    TVSS.μ = t -> (λ(t)*s+1 + sqrt((λ(t)*s+1)^2-4*s*(λ(t)-λ(t)*V)))/(2*s)
    TVSS.M = t -> QuadGK.quadgk(TVSS.μ, 0.0, t)[1]
    TVSS.M_interval = (x,y) -> QuadGK.quadgk(TVSS.μ, x, y)[1]
  elseif TVSS.control == "PD"
    TVSS.μ = t -> λ(t) + (V/s)
    TVSS.M = t -> TVAS.Λ(t)+t*(V/s)
    TVSS.M_interval = (x,y) -> TVAS.Λ_interval(x,y)+(y-x)*(V/s)
  elseif TVSS.control == "WS"
    TVSS.μ = t -> λ(t)+(λ(t)/2)*(sqrt(1+(4*V)/((s-1)*λ(t)))-1)
    TVSS.M = t -> QuadGK.quadgk(TVSS.μ, 0.0, t)[1]
    TVSS.M_interval = (x,y) -> QuadGK.quadgk(TVSS.μ, x, y)[1]
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
      push!(V, inverse_integrated_rate_function(TVAS.Λ, V[n], S, TVAS.table))
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

function next_event(system::TVGG1PS_queue, customer_pool::Array{Customer}, record::Record)
  if system.next_arrival_time == min(system.next_arrival_time, system.next_completion_time, system.next_virtual_completion_time, system.next_regular_recording)
    current_time = system.next_arrival_time
    # reduce workloads for each customer & virtual customer in system
    total_service_amount = 0.0
    if system.number_of_customers > 0 || system.number_of_virtual_customers > 0
      total_service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      individual_service_amount = total_service_amount/system.number_of_customers
      for customer in system.WIP
        customer.remaining_workload -= individual_service_amount # reduce workload of customers
      end
      individual_virtual_service_amount = total_service_amount/(system.number_of_customers+1)
      for virtual_customer in system.Virtual_WIP
        virtual_customer.remaining_workload -= individual_virtual_service_amount # reduce workload of virtual customers
      end
    end

    # insert a new customer in the system
    system.customer_arrival_counter += 1
    system.number_of_customers += 1
    push!(system.WIP, customer_pool[system.customer_arrival_counter])
    system.WIP[end].service_beginning_time = current_time
    # update simulational time
    system.sim_time = current_time

    # update next arrival time
    system.next_arrival_time = customer_pool[system.customer_arrival_counter+1].arrival_time

    # update next completion information
    if system.number_of_customers == 1
      system.next_completion_index = 1
      system.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, system.WIP[1].remaining_workload, system.TVSS.table)
    elseif system.number_of_customers > 1
      if system.WIP[system.next_completion_index].remaining_workload > system.WIP[end].remaining_workload
        system.next_completion_index = system.number_of_customers
      end
      nci = system.next_completion_index
      Q = system.number_of_customers
      system.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, Q*system.WIP[nci].remaining_workload, system.TVSS.table)
    end

    # update next virtual completion information
    if system.number_of_virtual_customers > 0
      nvci = system.next_virtual_completion_index
      Q = system.number_of_customers
      system.next_virtual_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, (Q+1)*system.Virtual_WIP[nvci].remaining_workload, system.TVSS.table)
    end

    total_workload = 0.0
    for customer in system.WIP
      total_workload += customer.remaining_workload
    end
  elseif system.next_completion_time == min(system.next_arrival_time, system.next_completion_time, system.next_virtual_completion_time, system.next_regular_recording)
    current_time = system.next_completion_time
    # reduce workloads for each customer & virtual customer in system
    total_service_amount = 0.0
    if system.number_of_customers > 1 || system.number_of_virtual_customers > 0 # if # of customer == 1, we can just remove the customer
      total_service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      individual_service_amount = total_service_amount/system.number_of_customers
      for customer in system.WIP
        customer.remaining_workload -= individual_service_amount # reduce workload of customers
      end
      individual_virtual_service_amount = total_service_amount/(system.number_of_customers+1)
      for virtual_customer in system.Virtual_WIP
        virtual_customer.remaining_workload -= individual_virtual_service_amount # reduce workload of virtual customers
      end
    end

    # save completed Customer time information
    nci = system.next_completion_index
    system.WIP[nci].completion_time = current_time
    system.WIP[nci].sojourn_time = current_time - system.WIP[nci].arrival_time
    system.WIP[nci].waiting_time = 0.0

    # copy the information of the completed customer to the customer pool
    customer_pool[system.WIP[nci].arrival_index] = system.WIP[nci]

    # remove the completing customer from the system
    deleteat!(system.WIP, nci)
    system.number_of_customers -= 1

    # update simulational time
    system.sim_time = current_time

    # update next completion information
    if system.number_of_customers == 0
      system.next_completion_index = 0
      system.next_completion_time = typemax(Float64)
    elseif system.number_of_customers == 1
      system.next_completion_index = 1
      system.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, system.WIP[1].remaining_workload, system.TVSS.table)
    elseif system.number_of_customers > 1
      shortest_remaining_workload = typemax(Float64)
      for i in 1:length(system.WIP)
        if shortest_remaining_workload > system.WIP[i].remaining_workload
          shortest_remaining_workload = system.WIP[i].remaining_workload
          system.next_completion_index = i
        end
      end
      nci = system.next_completion_index
      Q = system.number_of_customers
      system.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, Q*system.WIP[nci].remaining_workload, system.TVSS.table)
    end

    # update next virtual completion information
    if system.number_of_virtual_customers > 0
      nvci = system.next_virtual_completion_index
      Q = system.number_of_customers
      system.next_virtual_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, (Q+1)*system.Virtual_WIP[nvci].remaining_workload, system.TVSS.table)
    end
  elseif system.next_virtual_completion_time == min(system.next_arrival_time, system.next_completion_time, system.next_virtual_completion_time, system.next_regular_recording)
    current_time = system.next_virtual_completion_time
    # reduce workloads for each customer & virtual customer in system
    total_service_amount = 0.0
    if system.number_of_customers > 0 || system.number_of_virtual_customers > 1
      total_service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      individual_service_amount = total_service_amount/system.number_of_customers
      for customer in system.WIP
        customer.remaining_workload -= individual_service_amount # reduce workload of customers
      end
      individual_virtual_service_amount = total_service_amount/(system.number_of_customers+1)
      for virtual_customer in system.Virtual_WIP
        virtual_customer.remaining_workload -= individual_virtual_service_amount # reduce workload of virtual customers
      end
    end

    # Record the virtual sojourn time
    nvci = system.next_virtual_completion_index
    record.S[system.Virtual_WIP[nvci].time_index] = current_time - system.Virtual_WIP[nvci].arrival_time
    sojourn_time = record.S[system.Virtual_WIP[nvci].time_index]

    # remove the completing virtual customer from the system
    deleteat!(system.Virtual_WIP, system.next_virtual_completion_index)
    system.number_of_virtual_customers -= 1

    # update simulational time
    system.sim_time = current_time

    # update next virtual completion information
    if system.number_of_virtual_customers == 0
      system.next_virtual_completion_index = 0
      system.next_virtual_completion_time = typemax(Float64)
    elseif system.number_of_virtual_customers == 1
      system.next_virtual_completion_index = 1
      Q = system.number_of_customers
      system.next_virtual_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, (Q+1)*system.Virtual_WIP[1].remaining_workload, system.TVSS.table)
    elseif system.number_of_virtual_customers > 1
      shortest_remaining_workload = typemax(Float64)
      for i in 1:length(system.Virtual_WIP)
        if shortest_remaining_workload > system.Virtual_WIP[i].remaining_workload
          shortest_remaining_workload = system.Virtual_WIP[i].remaining_workload
          system.next_virtual_completion_index = i
        end
      end
      nvci = system.next_virtual_completion_index
      Q = system.number_of_customers
      system.next_virtual_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, (Q+1)*system.Virtual_WIP[nvci].remaining_workload, system.TVSS.table)
    end

  elseif system.next_regular_recording == min(system.next_arrival_time, system.next_completion_time, system.next_virtual_completion_time, system.next_regular_recording)
    current_time = system.next_regular_recording
    # reduce workloads for each customer & virtual customer in system
    total_service_amount = 0.0
    if system.number_of_customers > 0 || system.number_of_virtual_customers > 0
      total_service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      individual_service_amount = total_service_amount/system.number_of_customers
      for customer in system.WIP
        customer.remaining_workload -= individual_service_amount # reduce workload of customers
      end
      individual_virtual_service_amount = total_service_amount/(system.number_of_customers+1)
      for virtual_customer in system.Virtual_WIP
        virtual_customer.remaining_workload -= individual_virtual_service_amount # reduce workload of virtual customers
      end
    end

    # insert a new virtual customer in the system
    if system.TVSS.string_of_distribution != "H2"
      push!(system.Virtual_WIP, Virtual_Customer(system.time_index, rand(system.TVSS.workload_distribution), current_time, current_time, typemax(Float64), 0, 0))
    elseif system.TVSS.string_of_distribution == "H2"
      p1 = (5+sqrt(15))/10
      p2 = 1-p1
      θ1 = 1/(2*p1)
      θ2 = 1/(2*p2)
      push!(system.Virtual_WIP, Virtual_Customer(system.time_index, generate_Hyperexponential(p1,p2,θ1,θ2), current_time, current_time, typemax(Float64), 0, 0))
    end

    # increase # of virtual customer
    system.number_of_virtual_customers += 1

    # update simulational time & increase time index
    system.sim_time = current_time
    system.time_index += 1

    # update next regular recording time
    system.next_regular_recording += system.regular_recording_interval

    # update next virtual completion info
    if system.number_of_virtual_customers == 1
      system.next_virtual_completion_index = 1
      Q = system.number_of_customers
      system.next_virtual_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, (Q+1)*system.Virtual_WIP[1].remaining_workload, system.TVSS.table)
    else
      nvci = system.next_virtual_completion_index
      if system.Virtual_WIP[nvci].remaining_workload > system.Virtual_WIP[end].remaining_workload
        system.next_virtual_completion_index = system.number_of_virtual_customers
      end
      nvci = system.next_virtual_completion_index
      Q = system.number_of_customers
      system.next_virtual_completion_time = inverse_integrated_rate_function(system.TVSS.M, system.sim_time, (Q+1)*system.Virtual_WIP[nvci].remaining_workload, system.TVSS.table)
    end

    # record the number of customer, arrival process, temporal sojourn time
    push!(record.Q, system.number_of_customers)
    push!(record.A, system.customer_arrival_counter)
    push!(record.S, -1.0) # this will be re-recorded when the virtual customer is leaving
  end
end

function run_to_end(system::TVGG1PS_queue, customer_pool::Array{Customer}, record::Record, T::Float64, warm_up_time::Float64)
  while system.sim_time < T
    next_event(system, customer_pool, record)
  end
end

function do_experiment(queue::String, control::String, target::Float64, arrival::String, service::String, coeff::Tuple, T::Float64, N::Int64, record::Record)
  TVAS = Time_Varying_Arrival_Setting(coeff, arrival)
  TVSS = Time_Varying_Service_Setting(TVAS, control, target, service)
  set_distribution(TVAS)
  set_distribution(TVSS)
  set_service_rate_function(TVAS, TVSS)
  set_tables(TVAS,TVSS)
  file_num_in_queue = open("../../../logs/$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_queue_length.txt" , "w")
  file_virtual_sojourn_time = open("../../../logs/$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sojourn_time.txt" , "w")
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
    record.S = Float64[]
    record.A = Int64[]
    system = TVGG1PS_queue(TVAS, TVSS)
    customer_pool = generate_customer_pool(TVAS, TVSS, T*2.0)
    system.regular_recording_interval = T/1000
    system.next_arrival_time = customer_pool[1].arrival_time
    run_to_end(system, customer_pool, record, T*1.5, 0.0)
    writedlm(file_num_in_queue, transpose(record.Q)) # write record Q(t)
    writedlm(file_virtual_sojourn_time, transpose(record.S)) # write record W(t)
	end
  close(file_num_in_queue)
  close(file_virtual_sojourn_time)
end
