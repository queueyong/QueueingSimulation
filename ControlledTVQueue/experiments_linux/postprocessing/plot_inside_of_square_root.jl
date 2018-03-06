using PyPlot
figure(figsize = (15,5))
subplot(1,3,1)
title("s=0.1")
scatter(gmm.X[1,:],gmm.X[2,:]) # without clustering
subplot(1,3,2)
title("True Clusters")
for i in 1:length(mean)
 scatter(gmm.X[1,1+D*(i-1):D*i],gmm.X[2,1+D*(i-1):D*i])
end
Z = Points[]
for k in 1:gmm.K
 push!(Z,Points(Float64[],Float64[]))
end
for n in 1:gmm.N
 k = findmax(gmm.γ[n,:])[2]
 push!(Z[k].X, gmm.X[1,n])
 push!(Z[k].Y, gmm.X[2,n])
end
subplot(1,3,3)
title("Clustered Data Points")
for k in 1:gmm.K
 scatter(Z[k].X,Z[k].Y)
end
savefig("GMM clustering result.pdf")

V = 10.0
coeff = (1.0, 0.2, 0.001)
λ(t) = coeff[1]+coeff[2]*sin(coeff[3]*t)
f(t)=(λ(t)*s+1)^2-4*s*(λ(t)-λ(t)*V)
figure(figsize = (15,5))
subplot(1,3,1)
title("s=0.1, V=$V")
s = 0.1
Y = Float64[]
T = linspace(0.0,20*(1/coeff[3]),10000)
for t in T
 push!(Y,f(t))
end
plot(T,Y)
subplot(1,3,2)
title("s=1.0, V=$V")
s = 1.0
Y = Float64[]
T = linspace(0.0,20*(1/coeff[3]),10000)
for t in T
 push!(Y,f(t))
end
plot(T,Y)
subplot(1,3,3)
title("s=10.0, V=$V")
s = 10.0
Y = Float64[]
T = linspace(0.0,20*(1/coeff[3]),10000)
for t in T
 push!(Y,f(t))
end
plot(T,Y)
