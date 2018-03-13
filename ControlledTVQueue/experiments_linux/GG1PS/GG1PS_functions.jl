function scv(dist::Distribution)
  return var(dist)/mean(dist)^2
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

function string_to_dist(str::String) # All distributions have mean 1.0.
  if str == "ER1"
    return Erlang(2,9/20)
  elseif str == "ER2"
    return Erlang(2,9/10)
  elseif str == "LN1"
    return LogNormal(log(3*sqrt(3)/10),sqrt(log(3)))
  elseif str == "LN2"
    return LogNormal(log(3*sqrt(3)/5),sqrt(log(3)))
  elseif str == "ER3"
    return Erlang(1,1/2)
  elseif str == "ER4"
    return Erlang(2,1/2)
  elseif str == "LN3"
    return LogNormal(log(sqrt(3)/6), sqrt(log(3)))
  elseif str == "LN4"
    return LogNormal(log(sqrt(3)/3), sqrt(log(3)))
  elseif str == "ER0"
    return Erlang(2,1/2)
  elseif str == "LN0"
    return LogNormal(-(log(3)/2),sqrt(log(3)))
  end
end


function set_components(AS::Arrival_Setting, SS::Service_Setting)
  AS.interarrival_distribution = string_to_dist(AS.string_of_distribution)
  AS.λ = 1/mean(AS.interarrival_distribution)
  SS.workload_distribution = string_to_dist(SS.string_of_distribution)
end

function generate_NP(AS::Arrival_Setting, T::Float64)
    g_e = t -> 1-cdf(AS.interarrival_distribution,t)/mean(AS.interarrival_distribution)
    S = inverse_cdf(g_e , rand())
    V = [ S ]
    n = 1
    while V[n] < T
      push!(V, V[n] + rand(AS.interarrival_distribution))
      n += 1
    end
    return V
end

function generate_customer_pool(AS::Arrival_Setting, SS::Service_Setting, T::Float64)
    arrival_times = generate_NP(AS, T)
    Customers = Customer[]
    for i in 1:length(arrival_times)
      push!(Customers, Customer(i, rand(SS.workload_distribution), arrival_times[i], 0.0, typemax(Float64), 0.0, 0.0))
    end
    return Customers
end

function next_event(system::GG1PS_queue, customer_pool::Array{Customer}, record::Record)
  if system.next_arrival_time == min(system.next_arrival_time, system.next_completion_time, system.next_virtual_completion_time, system.next_regular_recording)
    current_time = system.next_arrival_time
    # reduce workloads for each customer & virtual customer in system
    total_service_amount = 0.0
    if system.number_of_customers > 0 || system.number_of_virtual_customers > 0
      total_service_amount = system.SS.μ*(current_time - system.sim_time)
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
      system.next_completion_time = system.sim_time + system.WIP[1].remaining_workload/system.SS.μ
    elseif system.number_of_customers > 1
      if system.WIP[system.next_completion_index].remaining_workload > system.WIP[end].remaining_workload
        system.next_completion_index = system.number_of_customers
      end
      nci = system.next_completion_index
      Q = system.number_of_customers
      system.next_completion_time = system.sim_time + Q*(system.WIP[nci].remaining_workload/system.SS.μ)
    end

    # update next virtual completion information
    if system.number_of_virtual_customers > 0
      nvci = system.next_virtual_completion_index
      Q = system.number_of_customers
      system.next_virtual_completion_time = system.sim_time + (Q+1)*(system.Virtual_WIP[nvci].remaining_workload/system.SS.μ)
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
      total_service_amount = system.SS.μ*(current_time - system.sim_time)
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

    # remove the completed customer from the system
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
      system.next_completion_time = system.sim_time + system.WIP[1].remaining_workload/system.SS.μ
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
      system.next_completion_time = system.sim_time + Q*(system.WIP[nci].remaining_workload/system.SS.μ)
    end

    # update next virtual completion information
    if system.number_of_virtual_customers > 0
      nvci = system.next_virtual_completion_index
      Q = system.number_of_customers
      system.next_virtual_completion_time = system.sim_time + (Q+1)*(system.Virtual_WIP[nvci].remaining_workload/system.SS.μ)
    end
  elseif system.next_virtual_completion_time == min(system.next_arrival_time, system.next_completion_time, system.next_virtual_completion_time, system.next_regular_recording)
    current_time = system.next_virtual_completion_time
    # reduce workloads for each customer & virtual customer in system
    total_service_amount = 0.0
    if system.number_of_customers > 0 || system.number_of_virtual_customers > 1
      total_service_amount = system.SS.μ*(current_time - system.sim_time)
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
      system.next_virtual_completion_time = system.sim_time + (Q+1)*(system.Virtual_WIP[1].remaining_workload/system.SS.μ)
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
      system.next_virtual_completion_time = system.sim_time + (Q+1)*(system.Virtual_WIP[nvci].remaining_workload/system.SS.μ)
    end

  elseif system.next_regular_recording == min(system.next_arrival_time, system.next_completion_time, system.next_virtual_completion_time, system.next_regular_recording)
    current_time = system.next_regular_recording
    # reduce workloads for each customer & virtual customer in system
    total_service_amount = 0.0
    if system.number_of_customers > 0 || system.number_of_virtual_customers > 0
      total_service_amount = system.SS.μ*(current_time - system.sim_time)
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
    push!(system.Virtual_WIP, Virtual_Customer(system.time_index, rand(system.SS.workload_distribution), current_time, current_time, typemax(Float64), 0, 0))

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
      system.next_virtual_completion_time = system.sim_time + (Q+1)*(system.Virtual_WIP[1].remaining_workload/system.SS.μ)
    else
      nvci = system.next_virtual_completion_index
      if system.Virtual_WIP[nvci].remaining_workload > system.Virtual_WIP[end].remaining_workload
        system.next_virtual_completion_index = system.number_of_virtual_customers
      end
      nvci = system.next_virtual_completion_index
      Q = system.number_of_customers
      system.next_virtual_completion_time = system.sim_time + (Q+1)*(system.Virtual_WIP[nvci].remaining_workload/system.SS.μ)
    end

    # record the number of customer, arrival process, temporal sojourn time
    push!(record.Q, system.number_of_customers)
    push!(record.A, system.customer_arrival_counter)
    push!(record.S, -1.0) # this will be re-recorded when the virtual customer is leaving
  end
end

function run_to_end(system::GG1PS_queue, customer_pool::Array{Customer}, record::Record, T::Float64, warm_up_time::Float64)
  while system.sim_time < T
    next_event(system, customer_pool, record)
  end
end

function do_experiment(queue::String, arrival::String, service::String, service_rate::Float64, T::Float64, N::Int64, record::Record)
  AS = Arrival_Setting(arrival)
  SS = Service_Setting(service, service_rate)
  set_components(AS, SS)
  file_queue_length = open("../logs/$(queue)/$(queue)_$(service_rate)_$(arrival)_$(service)_$(T)_$(N)_queue_length.txt" , "w")
  file_sojourn_time = open("../logs/$(queue)/$(queue)_$(service_rate)_$(arrival)_$(service)_$(T)_$(N)_sojourn_time.txt" , "w")
  regular_recording_interval = T/1000
  t = 0.0
  while t <= T
	   push!(record.T, t)
     t += regular_recording_interval
  end
  writedlm(file_queue_length, transpose(record.T)) # write time-axis
  writedlm(file_sojourn_time, transpose(record.T)) # write time-axis

  for n in 1:N
    println("Replication $n")
	  record.Q = Int64[]
    record.S = Float64[]
    record.A = Int64[]
    system = GG1PS_queue(AS, SS)
    customer_pool = generate_customer_pool(AS, SS, T*2.0)
    system.regular_recording_interval = T/1000
    system.next_arrival_time = customer_pool[1].arrival_time
    run_to_end(system, customer_pool, record, T*1.5, 0.0)
    writedlm(file_queue_length, transpose(record.Q)) # write record Q(t)
    writedlm(file_sojourn_time, transpose(record.S)) # write record W(t)
  end
  close(file_queue_length)
  close(file_sojourn_time)
  gc()
end
