using Distributions, PyPlot

function desired_scv_for_dist(str_dist::String)
 if str_dist == "ER"
  return 0.5
 elseif str_dist == "LN"
  return 2.0
 elseif str_dist == "EXP"
  return 1.0
 end
end


function SR(t, s, λ, β, arrival, service)
    Vps = (desired_scv_for_dist(arrival)+desired_scv_for_dist(service))/(1+desired_scv_for_dist(service))
    return  β*(λ(t) + (Vps/s))
end

function DM(t, s, λ, β, arrival, service)
    V = (desired_scv_for_dist(arrival)+desired_scv_for_dist(service))/2
    return ( β*(λ(t)*s+1) + sqrt( (β^2)*((λ(t)*s+1)^2) - (β^2)*4*s*(λ(t)-λ(t)*V) ) ) / (2*s)
end



# distributions
arrival = "ER"
service = "ER"

# arrival rate function
λ(t) = 1.0 + 0.2*sin(t)

# generate data
n = 100
t = linspace(0, 2*π, n)
s = linspace(0.5,10.0,n)

# make grids
xgrid = repmat(s',n,1)
ygrid = repmat(t,1,n)

# declare
z_SR = zeros(n,n)
z_DM = zeros(n,n)

for i in 1:n
    for j in 1:n
        z_SR[i,j] = SR(t[i],s[j],λ,1, arrival, service)
        z_DM[i,j] = DM(t[i],s[j],λ,1, arrival, service)
    end
end

#fig = figure("pyplot_surfaceplot",figsize=(10,10))
fig = figure()
PyPlot.PyObject(PyPlot.axes3D)
ax = fig[:add_subplot](1,1,1, projection = "3d")

#ax[:plot_surface](xgrid, ygrid, z_SR, rstride=2, edgecolors="k", cstride=2, color="red", alpha=0.7, linewidth=0.25, label=L"$\mu_{SR}(t;s)$")
#ax[:plot_surface](xgrid, ygrid, z_DM, rstride=2, edgecolors="k", cstride=2, color="blue", alpha=0.7, linewidth=0.25, label=L"$\mu_{DM}(t;s)$")
ax[:plot_surface](xgrid, ygrid, z_SR,  color="red", alpha=0.8, linewidth=1.5, label=L"$\mu_{SR}(t;s)$")
ax[:plot_surface](xgrid, ygrid, z_DM,  color="blue", alpha=0.8, linewidth=1.5, label=L"$\mu_{DM}(t;s)$")

#ax[:plot_surface](xgrid, ygrid, z_SR, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
xlabel("target response time (s)")
ylabel("time (t)")
legend(loc = "upper right", fontsize=12)
savefig("$(dirname(@__FILE__))/../plots/comparison of the two controls/comparison_$(arrival)_$(service).pdf")
