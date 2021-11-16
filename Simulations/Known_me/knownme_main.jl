using Pkg
Pkg.activate(".")

using Random, Distributions
using Plots
using LinearAlgebra, Optim
using Clustering

#Dir
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Simulations/Known_me/"
include(dirmain*"knownme_functions.jl")

pis=[2.0, 2.5, 3.0, 5.5];
p=[0.4,0.4,0.1, 0.1];

seed=10;
n=10000;

data=dgp(pis,p,n,seed);
histogram(data)
scatter(data,2.0*ones(n))
clusters=dbscan(data', 0.3)

func2(vars)=-loglike(vars,data)

basecluster, d, ntypes, pis_ini=f_basecluster(clusters,0)
param_ini=vcat(pis_ini,ones(length(pis_ini)-1)./length(pis_ini))
opt2 = optimize(func2, param_ini)
Optim.minimizer(opt2)
Optim.minimum(opt2)
func2(param_ini)

basecluster, d, ntypes, pis_ini=f_basecluster(clusters,1);
param_ini=vcat(pis_ini,ones(length(pis_ini)-1)./length(pis_ini));
opt2 = optimize(func2, param_ini);
Optim.minimizer(opt2)
Optim.minimum(opt2)

basecluster, d, ntypes, pis_ini=f_basecluster(clusters,15);
param_ini=vcat(pis_ini,ones(length(pis_ini)-1)./length(pis_ini));
opt2 = optimize(func2, param_ini);
Optim.minimizer(opt2)
Optim.minimum(opt2)
