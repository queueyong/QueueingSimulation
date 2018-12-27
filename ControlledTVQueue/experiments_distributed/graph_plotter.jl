using PyPlot
α = 1.0
β = 0.2
γ = 0.001
δ = 0.000025
T = 20000

fcn0 = t -> 1 + β*sin(γ*t)
fcn1 = t -> α + β*sin(γ*t)+β*sin(2*γ*t)
fcn1_string = "t -> α + β*sin(γ*t)+β*sin(2*γ*t)"
fcn2 = t -> α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t
fcn2_string = "t -> α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t"
fcn_dict = Dict("fcn1"=>fcn1_string, "fcn2"=>fcn2_string)
fcn_function_dict = Dict("fcn1"=>fcn1, "fcn2"=>fcn2)

x = linspace(0.0,T,10000)
y0 = fcn0(x)
y1 = fcn1(x)
y2 = [fcn2(t) for t in x]
plt = PyPlot
# Plotting arrival rates
plt.figure(figsize=(10,5))
plt.ylim(0,3.0)
plt.plot(x,y0,label="sin")
plt.plot(x,y1,label="fcn1")
plt.plot(x,y2,label="fcn2")
#plt.xlabel("Time (s)")
#plt.ylabel(L"$\lambda_i(s)$")
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc="upper right", fontsize=16)
plt.tight_layout()


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


function J(t::Float64, table::Array{Float64})
    index = Int64(floor(mod(t,(table[end-3]*table[end-2]))/table[end-1]))
    if index != 0
        return floor(t/(table[end-3]*table[end-2]))*table[end-2] + table[index]*table[end]
    else
        return floor(t/(table[end-3]*table[end-2]))*table[end-2]
    end
end

α = 1
β = 0.2
γ = 0.001
δ = 0.000025
THIS_FILE_PATH = dirname(@__FILE__)
queue = "TVGG1PS"
control = "SR"
target = 0.1
arrival = "ER"
service = "ER"
fcn = t -> α/2 + δ*t*sin(γ*t)+δ*t*sin(2*γ*t)+2*δ*t
fcn_type = "fcn2"
arrival_table = readdlm("$THIS_FILE_PATH/table/arrival/table_$(fcn_type)_arrival.txt")
service_table = readdlm("$THIS_FILE_PATH/table/service/$control/$target/$arrival,$service/table_$(fcn_type)_service_$(control)_$(arrival)_$(service).txt")


using Roots, QuadGK
T = 2000.0
α = 0.00025
β = 0.2
γ = 6*π/T
fcn0(t) = 1 + β*sin(γ*t)
fcn1(t) = 1 + β*sin(γ*t)+β*sin(2*γ*t)
fcn2(t) = α*t*sin(γ*t)+α*t*sin(2*γ*t)+2*α*t+0.5
fcn2(t) = 0.00025*t*sin(γ*t)+0.00025*t*sin(2*γ*t)+2*0.00025*t+0.5
table0 = table_for_any_periodic_J(fcn0, 2*π/γ)
table1 = table_for_any_periodic_J(fcn1, T/3)
table2 = table_for_any_periodic_J(fcn2, T)


t = 1242.0
x = J(t,arrival_table)
QuadGK.quadgk(fcn,0,x)[1]


t = 1242.0
x = J(t,table0)
QuadGK.quadgk(fcn0,0,x)[1]

t = 1242.0
x = J(t,table1)
QuadGK.quadgk(fcn1,0,x)[1]

t = 1242.0
x = J(t,table2)
QuadGK.quadgk(fcn2,0,x)[1]

using PyPlot
plt = PyPlot
arr = readdlm("$(dirname(@__FILE__))/TVGG1PS_0.1_PD_EXP_EXP_fnc1_2500.0_10000_sojourn_time.txt")
temp = Float64[]
for i in 1:size(arr)[2]
    push!(temp, mean(arr[2:end,i]))
end
temp
plt.figure()
plt.plot(arr[1,1:800],temp[1:800])

arr = readdlm("$(dirname(@__FILE__))/TVGG1PS_10.0_SR_LN_LN_fnc2_2000.0_10000_queue_length.txt")
mean(arr[:,1])
temp = Float64[]
for i in 1:size(arr)[2]
    push!(temp, mean(arr[:,i]))
end
plt.figure()
plt.plot(arr[1,:],temp)


#=
T = 20000

α = 0.00004
β = 0.2
γ = 6*π/T


fcn0(t) = 1 + β*sin(γ*t)
fcn1(t) = 1 + β*sin(γ*t)+β*sin(2*γ*t)
fcn2(t) = 0.00004*t*sin(γ*t)+0.00004*t*sin(2*γ*t)+2*0.00004*t+0.1
#fcn2(t) = 0.00004*t*sin(6*π/T*t)+0.00004*t*sin(2*6*π/T*t)+2*0.00004*t
x = linspace(0.0,T,10000)
y0 = fcn0(x)
y1 = fcn1(x)
y2 = [fcn2(t) for t in x]
plt = PyPlot
# Plotting arrival rates
plt.figure()
plt.ylim(0,2.5)

plt.plot(x,y0,label="sin")
plt.plot(x,y1,label="fcn1")
plt.plot(x,y2,label="fcn2")
#plt.xlabel("Time (s)")
#plt.ylabel(L"$\lambda_i(s)$")
plt.legend(loc="upper right")
plt.tight_layout()

=#
