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
seed=8;
n=10000;
data=dgp(pis,p,n,seed);
histogram(data)

#Solution
dens=epanechnikov
sup_length=2
T=5
npoints=100
τ=0.001
@time θ_hat, p_hat=mixture_dist(data, dens, sup_length, T, npoints, τ)

#Run code below to get simulations (can use pmap)
# function onesim(seed,pis,p,n,dens, sup_length, T, npoints, τ)
#     data=dgp(pis,p,n,seed);
#     θ_hat, p_hat=mixture_dist(data, dens, sup_length, T, npoints, τ);
#     return vcat(θ_hat,p_hat)
# end
#Result= map(seed->onesim(seed,pis,p,n,dens, sup_length, T, npoints, τ),1:1000);