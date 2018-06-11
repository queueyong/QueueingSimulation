include("./TVGG1PS_types.jl")

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

function set_tables(TVAS::Any, TVSS::Any)
    TVAS.table = table_for_sinusoidal_J((TVAS.α,TVAS.β,TVAS.γ))
    TVSS.table = table_for_any_periodic_J(TVSS.μ, 2*π/TVAS.γ)
end

function inverse_integrated_rate_function(Λ::Function, s::Float64, val::Float64, table::Array{Float64})
    return J(val+Λ(s),table)
end

function inverse_integrated_rate_function(Λ::Function, s::Float64, val::Float64)
    t = 0.0
    while Λ(t) < val + Λ(s)
        t += 0.00001
    end
    return t
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

function set_distribution(TVS::Any)
  if typeof(TVS) == Time_Varying_Arrival_Setting
    TVS.base_distribution = string_to_dist(TVS.string_of_distribution)
  elseif typeof(TVS) == Time_Varying_Service_Setting
    TVS.workload_distribution = string_to_dist(TVS.string_of_distribution)
  end
end

function set_service_rate_function(TVAS::Any, TVSS::Any)
  scv_arrival = scv(TVAS.base_distribution)
  scv_service = scv(TVSS.workload_distribution)
  β = mean(TVSS.workload_distribution)  # mean job size
  s = TVSS.target
  V = (scv_arrival+scv_service)/2
  Vps = (scv_arrival+scv_service)/(1+scv_service)
  λ = TVAS.λ
  if TVSS.control == "SR"
    TVSS.μ = t -> ( β*(λ(t)*s+1) + sqrt( (β^2)*((λ(t)*s+1)^2) - (β^2)*4*s*(λ(t)-λ(t)*V) ) ) / (2*s)
    TVSS.M = t -> QuadGK.quadgk(TVSS.μ, 0.0, t)[1]
    TVSS.M_interval = (x,y) -> QuadGK.quadgk(TVSS.μ, x, y)[1]
  elseif TVSS.control == "PD"
    TVSS.μ = t -> β*(λ(t) + (Vps/s))
    TVSS.M = t -> β*(TVAS.Λ(t)+t*(Vps/s))
    TVSS.M_interval = (x,y) -> β*(TVAS.Λ_interval(x,y)+(y-x)*(Vps/s))
  end
end

function generate_NHNP(TVAS::Time_Varying_Arrival_Setting, T::Float64)
    g_e = t -> 1-cdf(TVAS.base_distribution,t)/mean(TVAS.base_distribution)
    S = inverse_cdf(g_e , rand())
  	V = [ inverse_integrated_rate_function(TVAS.Λ, 0.0, S, TVAS.table) ]
    n = 1
    while V[n] < T
      push!(V, inverse_integrated_rate_function(TVAS.Λ, V[n], rand(TVAS.base_distribution), TVAS.table) )
      n += 1
    end
    return V
end

function generate_customer_pool(TVAS::Time_Varying_Arrival_Setting, TVSS::Time_Varying_Service_Setting, T::Float64)
    arrival_times = generate_NHNP(TVAS, T)
    Customers = Customer[]
    for i in 1:length(arrival_times)
      push!(Customers, Customer(i, rand(TVSS.workload_distribution), arrival_times[i], 0.0, typemax(Float64), 0.0, 0.0))
    end
    return Customers
end

function append_customers(customer_pool::Array{Customer}, TVAS::Time_Varying_Arrival_Setting, TVSS::Time_Varying_Service_Setting, T::Float64)
    former_length = length(customer_pool)
    former_last_arrival_time = customer_pool[end].arrival_time

    g_e = t -> 1-cdf(TVAS.base_distribution,t)/mean(TVAS.base_distribution)
    S = inverse_cdf(g_e , rand())
  	V = [ inverse_integrated_rate_function(TVAS.Λ, former_last_arrival_time, S, TVAS.table) ]
    n = 1
    while V[n] < former_last_arrival_time + T
      push!(V, inverse_integrated_rate_function(TVAS.Λ, V[n], rand(TVAS.base_distribution), TVAS.table) )
      n += 1
    end

    for n in 1:length(V)
        push!(customer_pool, Customer(former_length + n, rand(TVSS.workload_distribution), V[n], 0.0, typemax(Float64), 0.0, 0.0))
    end
end

# next_event for storing path (pre-simulation)
function next_event_pre(system::TVGG1PS_queue, customer_pool::Array{Customer}, record::Record)
  if system.next_arrival_time == min(system.next_arrival_time, system.next_completion_time, system.next_regular_recording)

    system.event = :A

    current_time = system.next_arrival_time
    # reduce workloads for each customer
    total_service_amount = 0.0
    if system.number_of_customers > 0
      total_service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      individual_service_amount = total_service_amount/system.number_of_customers
      for customer in system.WIP
        customer.remaining_workload -= individual_service_amount # reduce workload of customers
      end
    end

    # insert a new customer in the system
    system.customer_arrival_counter += 1
    system.number_of_customers += 1
    push!(system.WIP, deepcopy(customer_pool[system.customer_arrival_counter]))
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

    total_workload = 0.0
    for customer in system.WIP
      total_workload += customer.remaining_workload
    end
  elseif system.next_completion_time == min(system.next_arrival_time, system.next_completion_time, system.next_regular_recording)

    system.event = :C

    current_time = system.next_completion_time
    # reduce workloads for each customer & virtual customer in system
    total_service_amount = 0.0
    if system.number_of_customers == 1
        total_service_amount = system.TVSS.M_interval(system.sim_time, current_time)
        system.WIP[1].remaining_workload -= total_service_amount
    elseif system.number_of_customers > 1 # if # of customer == 1, we can just remove the customer
      total_service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      individual_service_amount = total_service_amount/system.number_of_customers
      for customer in system.WIP
        customer.remaining_workload -= individual_service_amount # reduce workload of customers
      end
    end

    # save completed Customer time information
    nci = system.next_completion_index
    system.WIP[nci].completion_time = current_time
    system.WIP[nci].sojourn_time = current_time - system.WIP[nci].arrival_time
    system.WIP[nci].waiting_time = 0.0

    # copy the information of the completed customer to the customer pool
    #customer_pool[system.WIP[nci].arrival_index] = deepcopy(system.WIP[nci])

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
  elseif system.next_regular_recording == min(system.next_arrival_time, system.next_completion_time, system.next_regular_recording)

    system.event = :R

    current_time = system.next_regular_recording
    # reduce workloads for each customer in system
    total_service_amount = 0.0
    if system.number_of_customers > 0
      total_service_amount = system.TVSS.M_interval(system.sim_time, current_time)
      individual_service_amount = total_service_amount/system.number_of_customers
      for customer in system.WIP
        customer.remaining_workload -= individual_service_amount # reduce workload of customers
      end
    end

    # update simulational time & increase time index
    system.sim_time = current_time
    system.time_index += 1

    # update next regular recording time
    system.next_regular_recording += system.regular_recording_interval

    # store current system information (Customer array)
    #TODO make sure that the copy is well done (deep, shallow issue)
    record.P[current_time] = TVGG1PS_queue_path(system)
    push!(record.Q, system.number_of_customers)
  end
end

function next_event_vc(system::TVGG1PS_queue, path::TVGG1PS_queue_path, customer_pool::Array{Customer})
    if path.next_arrival_time == min(path.next_arrival_time, path.next_completion_time)
      path.event = :A

      current_time = path.next_arrival_time
      # reduce workloads for each customer
      total_service_amount = 0.0
      if path.number_of_customers > 0
        total_service_amount = system.TVSS.M_interval(path.sim_time, current_time)
        individual_service_amount = total_service_amount/path.number_of_customers
        for customer in path.WIP
          customer.remaining_workload -= individual_service_amount # reduce workload of customers
        end
      end

      # insert a new customer in the system
      path.customer_arrival_counter += 1
      path.number_of_customers += 1
      push!(path.WIP, deepcopy(customer_pool[path.customer_arrival_counter]))
      path.WIP[end].service_beginning_time = current_time

      # update simulational time
      path.sim_time = current_time

      # check if customer_pool is insufficient & generate new customers if so
      if path.customer_arrival_counter == path.size_customer_pool
          println("appending customers")
          append_customers(customer_pool, system.TVAS, system.TVSS, 1000.0)
          path.size_customer_pool = length(customer_pool)
      end

      # update next arrival time
      path.next_arrival_time = customer_pool[path.customer_arrival_counter+1].arrival_time

      # update next completion information
      if path.number_of_customers == 1
        path.next_completion_index = 1
        path.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, path.sim_time, path.WIP[1].remaining_workload, system.TVSS.table)
      elseif path.number_of_customers > 1
        if path.WIP[path.next_completion_index].remaining_workload > path.WIP[end].remaining_workload
          path.next_completion_index = path.number_of_customers
        end
        nci = path.next_completion_index
        Q = path.number_of_customers
        path.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, path.sim_time, Q*path.WIP[nci].remaining_workload, system.TVSS.table)
      end

      total_workload = 0.0
      for customer in path.WIP
        total_workload += customer.remaining_workload
      end
    elseif path.next_completion_time == min(path.next_arrival_time, path.next_completion_time)
      path.event = :C

      current_time = path.next_completion_time
      # reduce workloads for each customer & virtual customer in system
      total_service_amount = 0.0
      if path.number_of_customers == 1
          total_service_amount = system.TVSS.M_interval(path.sim_time, current_time)
          path.WIP[1].remaining_workload -= total_service_amount
      elseif path.number_of_customers > 1 # if # of customer == 1, we can just remove the customer
        total_service_amount = system.TVSS.M_interval(path.sim_time, current_time)
        individual_service_amount = total_service_amount/path.number_of_customers
        for customer in path.WIP
          customer.remaining_workload -= individual_service_amount # reduce workload of customers
        end
      end

      # save completed Customer time information
      nci = path.next_completion_index
      path.WIP[nci].completion_time = current_time
      path.WIP[nci].sojourn_time = current_time - path.WIP[nci].arrival_time
      path.WIP[nci].waiting_time = 0.0

      # remove the completing customer from the system
      deleteat!(path.WIP, nci)
      path.number_of_customers -= 1

      # update simulational time
      path.sim_time = current_time

      # update next completion information
      if path.number_of_customers == 0
        path.next_completion_index = 0
        path.next_completion_time = typemax(Float64)
      elseif path.number_of_customers == 1
        path.next_completion_index = 1
        path.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, path.sim_time, path.WIP[1].remaining_workload, system.TVSS.table)
      elseif path.number_of_customers > 1
        shortest_remaining_workload = typemax(Float64)
        for i in 1:length(path.WIP)
          if shortest_remaining_workload > path.WIP[i].remaining_workload
            shortest_remaining_workload = path.WIP[i].remaining_workload
            path.next_completion_index = i
          end
        end
        nci = path.next_completion_index
        Q = path.number_of_customers
        path.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, path.sim_time, Q*path.WIP[nci].remaining_workload, system.TVSS.table)
      end
    end
end

function storePath(system::TVGG1PS_queue, customer_pool::Array{Customer}, record::Record, T::Float64)
    while system.sim_time < T
      next_event_pre(system, customer_pool, record)
    end
#    return record.P
end

function run_to_end(system::TVGG1PS_queue, customer_pool::Array{Customer}, record::Record)
    for t in record.T
        push!(record.S, simulateSojournTime(system, customer_pool, record.P[t]))
    end
end

## run simulation proximally from recording epoch t to simulate/record virtual sojourn time at t
function simulateSojournTime(system::TVGG1PS_queue, customer_pool::Array{Customer}, path::TVGG1PS_queue_path)
#    println("simlating at $(path.sim_time)")

    # create and insert virtual customer
    vc = Customer(-1, rand(system.TVSS.workload_distribution), path.sim_time, path.sim_time, typemax(Float64), -1.0, 0.0)
    push!(path.WIP, vc)

    # update system status
    if path.sim_time != 0.0
        ## increase number of customer
        path.number_of_customers += 1
        ## update next completion information
        if path.number_of_customers == 1
          path.next_completion_index = 1
          path.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, path.sim_time, path.WIP[1].remaining_workload, system.TVSS.table)
        elseif path.number_of_customers > 1
          if path.WIP[path.next_completion_index].remaining_workload > path.WIP[end].remaining_workload
            path.next_completion_index = path.number_of_customers
          end
          nci = path.next_completion_index
          Q = path.number_of_customers
          path.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, path.sim_time, Q*path.WIP[nci].remaining_workload, system.TVSS.table)
        end
    else
        ## increase number of customer
        path.number_of_customers += 1
        ## update next completion information
        path.next_completion_index = 1
        path.next_completion_time = inverse_integrated_rate_function(system.TVSS.M, path.sim_time, path.WIP[1].remaining_workload)
    end

    # update size of customer pool
    path.size_customer_pool = length(customer_pool)

    # run simulation till the virtual customer leaves
    while vc.sojourn_time == -1.0
        next_event_vc(system, path, customer_pool)
    end

    return vc.sojourn_time
end

function defaultize(system::TVGG1PS_queue)
    system.WIP = Customer[]
    system.number_of_customers = 0
    system.sim_time = 0.0
    system.next_arrival_time = typemax(Float64)
    system.next_completion_time = typemax(Float64)
    system.next_regular_recording = 0.0
    system.next_completion_index = 0
    system.time_index = 1
    system.customer_arrival_counter = 0
    system.event = :A
    system.size_customer_pool = 0
end

## To fully simulate a queueing system (records Q(t), S(t))
function do_experiment(queue::String, control::String, target::Float64, arrival::String, service::String, coeff::Tuple, T::Float64, N::Int64, record::Record, LOG_PATH::String)
    TVAS = Time_Varying_Arrival_Setting(coeff, arrival)
    TVSS = Time_Varying_Service_Setting(TVAS, control, target, service)
    set_distribution(TVAS)
    set_distribution(TVSS)
    set_service_rate_function(TVAS, TVSS)
    set_tables(TVAS,TVSS)
    file_num_in_queue = open("$LOG_PATH/$(queue)/$(target)/$(arrival),$(service)/$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_queue_length.txt" , "w")
    file_virtual_sojourn_time = open("$LOG_PATH/$(queue)/$(target)/$(arrival),$(service)/$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N)_sojourn_time.txt" , "w")
    running_mark = open("$LOG_PATH/$(queue)/$(target)/$(arrival),$(service)/running_$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N).txt", "w")

    system = TVGG1PS_queue(TVAS, TVSS)
    system.regular_recording_interval = T/1000
    record.T = [t for t in 0.0:system.regular_recording_interval:T]
    writedlm(file_num_in_queue, transpose(record.T)) # write time-axis
    writedlm(file_virtual_sojourn_time, transpose(record.T)) # write time-axis

    for n in 1:N
      write(running_mark, "Replication $n")
      println("Replication $n")
      record.Q = Int64[]
      record.S = Float64[]
      record.P = Dict{Float64, TVGG1PS_queue_path}()
      customer_pool = generate_customer_pool(TVAS, TVSS, T*1.2)
      defaultize(system)
      system.next_arrival_time = customer_pool[1].arrival_time
      system.size_customer_pool = length(customer_pool)
      storePath(system, customer_pool, record, T)
      run_to_end(system, customer_pool, record)
      writedlm(file_num_in_queue, transpose(record.Q)) # write record Q(t)
      writedlm(file_virtual_sojourn_time, transpose(record.S)) # write record W(t)
    end
    close(running_mark)
    close(file_num_in_queue)
    close(file_virtual_sojourn_time)

    rm("$LOG_PATH/$(queue)/$(target)/$(arrival),$(service)/running_$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N).txt")
    completion_mark = open("$LOG_PATH/$(queue)/$(target)/$(arrival),$(service)/completed_$(queue)_$(target)_$(control)_$(arrival)_$(service)_$(coeff[3])_$(T)_$(N).txt", "w")
    close(completion_mark)

    gc()
end

#=

# revised simulation functions
queue = "TVGG1PS"
control = "SR"
target = 10.0
arrival = "EXP"
service = "EXP"
coeff = (1.0, 0.2, 0.1)
T = 2000.0
N = 10000
record = Record()


TVAS = Time_Varying_Arrival_Setting(coeff, arrival)
TVSS = Time_Varying_Service_Setting(TVAS, control, target, service)
set_distribution(TVAS)
set_distribution(TVSS)
set_service_rate_function(TVAS, TVSS)
set_tables(TVAS,TVSS)



system = TVGG1PS_queue(TVAS, TVSS)
system.regular_recording_interval = T/1000
record.T = [t for t in 0.0:system.regular_recording_interval:T]

record.Q = Int64[]
record.S = Float64[]
record.P = Dict{Float64, TVGG1PS_queue_path}()
customer_pool = generate_customer_pool(TVAS, TVSS, T*1.2)
defaultize(system)
system.next_arrival_time = customer_pool[1].arrival_time
system.size_customer_pool = length(customer_pool)
storePath(system, customer_pool, record, T)
run_to_end(system, customer_pool, record)


writedlm(file_num_in_queue, transpose(record.Q)) # write record Q(t)
writedlm(file_virtual_sojourn_time, transpose(record.S)) # write recor


system = TVGG1PS_queue(TVAS, TVSS)
system.regular_recording_interval = T/1000
record = Record()
record.T = [t for t in 0.0:system.regular_recording_interval:T]
customer_pool = generate_customer_pool(system.TVAS, system.TVSS, T*1.5)
system.next_arrival_time = customer_pool[1].arrival_time
system.size_customer_pool = length(customer_pool)
path = storePath(system, customer_pool, record, T)
run_to_end(system, customer_pool, record, path)

record.S

do_experiment(queue, control, target, arrival, service, coeff, T, N, record)

@enter do_experiment(queue, control, target, arrival, service, coeff, T, N, record)

=#
