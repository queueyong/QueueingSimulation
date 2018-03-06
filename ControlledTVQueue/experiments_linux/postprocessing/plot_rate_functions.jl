cd(dirname(Base.source_path()))
using PyPlot

queue = "TVGG1PS"
coeff = (1.0, 0.2, 0.001)
control = "SR"
arrival, service = "E2", "E2"
scv_arrival = 0.0
scv_service = 0.0
if arrival == "EXP"
  scv_arrival = 1.0
elseif arrival == "H2"
  scv_arrival = 4.0
elseif arrival == "E2"
  scv_arrival = 0.5
end
if service == "EXP"
  scv_service = 1.0
elseif service == "H2"
  scv_service = 4.0
elseif service == "E2"
  scv_service = 0.5
end

s = 1.0
V = 0.0 #temporary
if queue == "TVGG1"
  V = (scv_arrival+scv_service)/2
else
  V = (scv_arrival +scv_service)/(1+scv_service)
end

λ(t) = coeff[1]+coeff[2]*sin(coeff[3]*t)
μ_SR(t) = (λ(t)*s+1 + sqrt((λ(t)*s+1)^2-4*s*(λ(t)-λ(t)*V)))/(2*s)
μ_PD(t) = λ(t) + (V/s)

plt = PyPlot
plt.figure()
plt.title("λ(t)=$(coeff[1])+$(coeff[2])sin($(coeff[3])t), s=$s, arrival:$arrival/service:$service (V=$V)")
X = linspace(0.0,20*(1/coeff[3]),10000)
λ_array = Float64[]
μ_SR_array = Float64[]
μ_PD_array = Float64[]
for x in X
  push!(λ_array, λ(x))
  push!(μ_SR_array, μ_SR(x))
  push!(μ_PD_array, μ_PD(x))
end
plt.plot(X,λ_array,label="λ(t)")
plt.plot(X,μ_SR_array,label=L"$μ_{SR}(t)$")
plt.plot(X,μ_PD_array,label=L"$μ_{PD}(t)$")
plt.xlim(0.0,20*(1/coeff[3]),10000)
plt.ylim(0.0,1.2*max(maximum(λ_array),maximum(μ_SR_array),maximum(μ_PD_array)))
plt.legend()
plt.savefig("../plots/rate_functions_$(queue)_γ_$(coeff[3])_s_$(s)_$(arrival)_$(service).pdf")

###############################################################################################33
coeff = (1.0, 0.2, 0.001)
s = 0.1
V = 4.0
λ(t) = coeff[1]+coeff[2]*sin(coeff[3]*t)
μ_plus(t) = (λ(t)*s+1 + sqrt((λ(t)*s+1)^2-4*s*(λ(t)-λ(t)*V)))/(2*s)
μ_minus(t) = (λ(t)*s+1 - sqrt((λ(t)*s+1)^2-4*s*(λ(t)-λ(t)*V)))/(2*s)
μ_SR(t) = λ(t)/2-(1/(2*s))+sqrt(((s*λ(t)+1)^2)-4*s*(1-V)*λ(t))/(2*s)
plt = PyPlot
plt.figure()
X = linspace(0.0,20*(1/coeff[3]),10000)
λ_array = Float64[]
μ_plus_array = Float64[]
μ_minus_array = Float64[]
μ_SR_array = Float64[]
for x in X
  push!(λ_array, λ(x))
  push!(μ_plus_array, μ_plus(x))
  push!(μ_minus_array, μ_minus(x))
  push!(μ_SR_array, μ_SR(x))
end
plt.plot(X,λ_array,label="λ(t)")
plt.plot(X,μ_plus_array,label=L"$μ_{+}(t)$")
plt.plot(X,μ_minus_array,label=L"$μ_{-}(t)$")
plt.plot(X,μ_SR_array,label=L"$μ_{SR}(t)$")
plt.xlim(0.0,20*(1/coeff[3]),10000)
#plt.ylim(0.0,1.2*max(maximum(λ_array),maximum(μ_plus_array),maximum(μ_minus_array)))
plt.legend()
plt.savefig("../plots/mu_plus_minus.pdf")
