using Roots, QuadGK
include("$(dirname(@__FILE__))/../functions_distribution.jl")

function service_rate_function(control::String, arrival_function::Function, target::Float64, arrival::String, service::String)
  dist_arrival = string_to_dist(arrival)
  dist_service = string_to_dist(service)
  scv_arrival = scv(dist_arrival)
  scv_service = scv(dist_service)
  β = mean(dist_service)  # mean job size
  s = target
  V = (scv_arrival+scv_service)/2
  Vps = (scv_arrival+scv_service)/(1+scv_service)
  λ = arrival_function
  if control == "SR"
    return t -> ( β*(λ(t)*s+1) + sqrt( (β^2)*((λ(t)*s+1)^2) - (β^2)*4*s*(λ(t)-λ(t)*V) ) ) / (2*s)
  elseif control == "PD"
    return t -> β*(λ(t) + (Vps/s))
  end
end

function table_for_any_periodic_J(f::Function, T::Float64)  # standard_J: λ(t)=1+βsin(t)
    temp = Float64[]
    X = 0.0:T/1000:T
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


T = 20000.0
α = 1
β = 0.2
γ = 0.001
δ = 0.000025

fcn0(t) = α + β*sin(γ*t)
fcn1(t) = α + β*sin(γ*t)+β*sin(2*γ*t)
fcn2(t) = α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t

for fcn in [fcn0, fcn1, fcn2]
    # generate arrival table
    table_period = 0.0
    if fcn != fcn2
        table_period = 2*π/γ
    else
        table_period = 22000.0
    end
    table = table_for_any_periodic_J(fcn, table_period)
    file = open("$(dirname(@__FILE__))/arrival/table_$(fcn)_arrival.txt","w")
    writedlm(file, table)
    close(file)

    # generate service table
    for control in ["PD", "SR"]
        for target in [0.1, 10.0]
            for (arrival,service) in [("EXP","EXP"), ("ER","ER"), ("LN","LN"), ("ER","LN"), ("LN","ER")]
                μ = service_rate_function(control, fcn, target, arrival, service)
                table = table_for_any_periodic_J(μ, table_period)
                file = open("$(dirname(@__FILE__))/service/$control/$target/$arrival,$service/table_$(fcn)_service_$(control)_$(arrival)_$(service).txt","w")
                writedlm(file, table)
                close(file)
            end
        end
    end
end
