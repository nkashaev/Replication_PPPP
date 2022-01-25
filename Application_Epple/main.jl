using Pkg
Pkg.activate(".")

using LinearAlgebra
using Distributions, Statistics, Random
using Clustering
using LowRankApprox, Printf
using JuMP, KNITRO
using KDEstimation
using CSV, DataFrames
#################################### Dir ###############################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Application_Epple/"
dirdata=dirmain*"data"
dirresults=dirmain*"results"
#################################### Functions ###############################
include(dirmain*"application_functions.jl")
include(dirmain*"mixSQP.jl")

#################################### Data #################################### 
nc=10
data=CSV.read(dirdata*"/data_cleaned_$(nc).csv", DataFrame)
pl=data[:,1]; v=data[:,2]; market=data[:,3]; zip=data[:,4];
#################################### Deconvolution of mismeasured values #################################### 
### Estimation of the pdf of the measurement error
if nc==10
    #Markets 8 and 9 have very similar upper clusters
    #We will demean them, pull them together and estimate the measurement error density
    dataforme1=sort(v[market.==8])[379:end]
    dataforme2=sort(v[market.==9])[211:end]
else
    #Markets 6 and 7 have very similar upper clusters
    #We will demean them, pull them together and estimate the measurement error density
    dataforme1=sort(v[market.==6])[465:end]
    dataforme2=sort(v[market.==7])[241:end]
end

dataforme=sort(vcat(dataforme1 .-mean(dataforme1),dataforme2 .-mean(dataforme2)))
# Kernel Estimation
kerfunB=Epanechnikov
kerfunK=epanechnikov
dens, sup_length=mixed_density_app(dataforme, kerfunB,kerfunK)
# rn=-10:0.1:10
# plot(rn,dens.(rn))
### Estimation of v(pl,e)
if nc==10
    T=5
else
    T=4 # Number of typed
end
npoints=100 # Number of discritizations on the 1st step
τ=0.0001 #Threshold for trimming
twostep=true #Set to true to use 
Theta=Array{Vector{Float64}}(undef,nc)
P=Array{Vector{Float64}}(undef,nc)
model = Model(with_optimizer(KNITRO.Optimizer))
for m in 1:nc
    data=sort(v[market.==m])
    θ_hat, p_hat=mixture_dist1(data, dens, sup_length, T, npoints, τ)
    if twostep
        θ_hat, p_hat=mixture_dist2(data, dens, vcat(θ_hat, p_hat), model)
    end
    Theta[m]=θ_hat
    P[m]=p_hat
end

#################################### Estimation of proxy function #################################### 
##Computing average v and pl per market
vbar=zeros(nc); plbar=zeros(nc);
for m in 1:nc
    vbar[m]=mean(v[market.==m]); plbar[m]=mean(pl[market.==m]);
end

f_po=g_proxy(vbar,plbar,3);

#################################### Estimation of the suply function #################################### 
# Theta=Theta[3:end] #First market has only 4 types
# vbar=vbar[3:end]
yo=zeros(T,length(Theta));
po=zeros(length(Theta));

for m in 1:length(Theta)
    po[m]=f_po(vbar[m])
    yo[:,m]=sort(Theta[m])./po[m]
end
plot(yo')
# [(log.(po[m]),log.(yo[4,m]')) for m=2:8]
#Saving results
CSV.write(dirresults*"/output_level_$(nc)_$(T).csv", DataFrame(yo,:auto))
CSV.write(dirresults*"/output_price_$(nc)_$(T).csv", DataFrame(po',:auto))
