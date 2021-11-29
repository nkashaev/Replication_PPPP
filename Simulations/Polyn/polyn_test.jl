using Pkg
Pkg.activate(".")

using Random, Distributions, LinearAlgebra
using JuMP, SCS
using Plots, Clustering
using LowRankApprox, Printf
using Optim
using KDEstimation
#Dir
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Simulations/Polyn/"
#Functions
include(dirmain*"polyn_functions.jl")
include(dirmain*"mixSQP.jl")

# Parameter of DGP
pis=[2.0, 2.5, 3.0, 3.5, 6.0];
p=[0.15,0.25,0.1, 0.2, 0.3];
seed=8;
n=1000;
data=dgp(pis,p,n,seed);
histogram(data)

#Kernel estimation
kerfun=Epanechnikov
# kerfun=Normal
# kerfun=Logistic
# kerfun=SymTriangularDist
cluster_distance=0.3
dens, sup_length=mixed_density(data, kerfun, cluster_distance)
plot(-1.5:0.001:1.5,dens.(-1.5:0.001:1.5))
#Solution
T=5
npoints=100
τ=0.001
twostep=false
@time θ_hat, p_hat=mixture_dist(data, dens, sup_length, T, npoints, τ, twostep)

#Run code below to get simulations (can use pmap)
# function onesim(seed,pis,p,n,dens, sup_length, T, npoints, τ)
#     data=dgp(pis,p,n,seed);
#     θ_hat, p_hat=mixture_dist(data, dens, sup_length, T, npoints, τ);
#     return vcat(θ_hat,p_hat)
# end
#Result= map(seed->onesim(seed,pis,p,n,dens, sup_length, T, npoints, τ),1:1000);