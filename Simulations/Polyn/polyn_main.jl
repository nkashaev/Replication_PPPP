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
n=100;

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



using JuMP
using KNITRO
using GLPK

theta=[1.1+k*0.05 for k in 0:100]
A=zeros(n,length(theta))
for i in 1:n, j in 1:length(theta)
    A[i,j]=epanechnikov(data[i] - theta[j])
end

model = Model(optimizer_with_attributes(ECOS.Optimizer, "printlevel" => 0))
@variable(model, v[1:n]);
@variable(model, t[1:n]>= 0.0);
Obj=@expression(model,sum(v[i] for i in 1:n));
@objective(model, Max, Obj)
@constraint(model, c1[j=1:length(theta)], sum(A[i,j]*t[i] for i=1:n) <= n);
@constraint(model, con[i = 1:n], [v[i], 1.0, t[i]] in MOI.ExponentialCone())
print(model)
optimize!(model)

ss