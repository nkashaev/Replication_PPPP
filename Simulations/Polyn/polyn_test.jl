using Pkg
Pkg.activate(".")

using Random, Distributions, LinearAlgebra
using JuMP, SCS
using Plots, Clustering
using LowRankApprox, Printf
using Optim
#Dir
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Simulations/Polyn/"
#Functions
include(dirmain*"polyn_functions.jl")
include(dirmain*"mixSQP.jl")

# Parameter of DGP
pis=[2.0, 2.5, 3.0, 3.5, 5.5];
p=[0.25,0.35,0.1, 0.2, 0.1];
seed=10;
n=300;
data=dgp(pis,p,n,seed);
histogram(data)


## SQP method
npoints=100
θ=[minimum(data)+k*(maximum(data)-minimum(data))/npoints for k in 2:npoints-2]
out = mixSQP(MatA(θ,data),maxiter = 10,verbose = false);
out1 = mixSQP(MatA(θ,data),x = ones(length(θ))/length(θ));

#Taking solutions with large enough weigths
τ=0.01
p_raw=out1["x"]
θ_raw=θ[p_raw .>τ]
p_raw=p_raw[p_raw .>τ]

#Clustering θ to T types
T=5
cl_θ=kmeans(θ_raw', T; weights=p_raw)
θ_hat=sort(cl_θ.centers[:]) #Intial guess for θ
out2 = mixSQP(MatA(θ_hat,data),x = ones(length(θ_hat))/length(θ_hat));
p_hat=out2["x"] #Intial guess for p
#Support adjustment
θ_hat=supportadj!(θ_hat,p_hat,data) 
#Final optimization
param_ini=vcat(θ_hat,log.(ones(length(θ_hat))[1:end-1]/length(θ_hat))) # Log is taken since we use exp(variable) to guarantee that p>0
opt2 = optimize(vars->-loglike(vars,data), param_ini, iterations=10^5)
sol=Optim.minimizer(opt2)
#Checking if the true parameter value gives a better objective function
println("The solution is better than the true parameter value --",(-loglike(vcat(pis,log.(p[1:end-1])),data)>Optim.minimum(opt2)))
#Final estimates
θ_hat=sol[1:length(θ_hat)]
p_hat=vcat(exp.(sol[length(θ_hat)+1:end]),1.0-sum(exp.(sol[length(θ_hat)+1:end])))