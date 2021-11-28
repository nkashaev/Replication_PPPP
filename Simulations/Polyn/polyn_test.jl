using Pkg
Pkg.activate(".")

using Random, Distributions, LinearAlgebra
using JuMP, SCS
using Plots, Clustering
using LowRankApprox, Printf
#Dir
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Simulations/Polyn/"
#Functions
include(dirmain*"polyn_functions.jl")
include(dirmain*"mixSQP.jl")

# Parameter of DGP
pis=[2.0, 2.5, 3.0, 5.5];
p=[0.4,0.4,0.1,0.1];
seed=10;
n=100000;
data=dgp(pis,p,n,seed);
#histogram(data)


## SQP method
npoints=100
θ=[minimum(data)+k*(maximum(data)-minimum(data))/dtheta for k in 2:dtheta-2]
out = mixSQP(MatA(θ,data),maxiter = 10,verbose = false);
out1 = mixSQP(MatA(θ,data),x = ones(length(θ))/length(θ));

#Taking solutions with large enough weigths
τ=0.01
p_raw=out1["x"]
θ_raw=θ[p_raw .>τ]
p_raw=p_raw[p_raw .>τ]

#Clustering θ to T types
T=4
cl_θ=kmeans(θ_raw', T; weights=p_hat_raw)
θ_hat=sort(cl_θ.centers[:]) #Intial guess for θ
out2 = mixSQP(MatA(θ_hat,data),x = ones(length(θ_hat))/length(θ_hat));
p_hat=out["x"] #Intial guess for p

using Optim
func2(vars)=-loglike(vars,data)
#∇func2!(G,vars)=-∇loglike!(G,vars,data)
param_ini=vcat(θ_hat,log.([0.2,0.2,0.2]))
param_ini[1]=minimum(data)+0.99
func2(param_ini)

opt2 = optimize(func2, param_ini)
sol=Optim.minimizer(opt2)
val=Optim.minimum(opt2)

#func2(vcat(pis,log.(p[1:end-1])))
#Final estimates
θ_hat=sol[1:length(θ_hat)]
p_hat=vcat(exp.(sol[length(θ_hat)+1:end]),1.0-sum(exp.(sol[length(θ_hat)+1:end])))