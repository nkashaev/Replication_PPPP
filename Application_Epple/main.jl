using Pkg
Pkg.activate(".")

using LinearAlgebra, DelimitedFiles
using Distributions, Statistics, Random
using Plots
using Clustering#, ParallelKMeans
using LowRankApprox, Printf
using JuMP, KNITRO, Optim
using KDEstimation
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
EppleData=readdlm(dirdata*"/Pittsburgh_post1995.txt", ',', Float64, '\n',header=true);
Data=EppleData[1];
zipcodedata_raw=readdlm(dirdata*"/US Zip Codes from 2013 Government Data.txt", ',', Float64, '\n',header=true)[1];
#Removing location, pl, and v outliers. Also removing houses with high pl and value 
Noutliersxy= (-80.3.<Data[:,9].<-79.67).*(40.25.<Data[:,10]);
Noutliersplv= (Data[:,2].<7).*(Data[:,3].<55);
println("The percent of dropped observations is ", 100*(1.0-mean(Noutliersplv))) # percent of dropped observations

x=Data[Noutliersxy.*Noutliersplv,9]; y=Data[Noutliersxy.*Noutliersplv,10]; #coordinates
pl=Data[Noutliersxy.*Noutliersplv,2]; v=Data[Noutliersxy.*Noutliersplv,3]; #price of land and value
N=length(pl);
#Zipcodes
Z=(minimum(y).<=zipcodedata_raw[:,2].<=maximum(y)).*(minimum(x).<=zipcodedata_raw[:,3].<=maximum(x))
zip_coord=zipcodedata_raw[Z,:]
zip=zeros(N)
for i in 1:N
    D=(x[i] .- zip_coord[:,3]).^2 .+ (y[i] .- zip_coord[:,2]).^2
    zip[i]=zip_coord[argmin(D),1]
end
#scatter(zip, pl, color=:lightrainbow, legend=false)
scatter(pl,v)

#################################### Preparing Markets #################################### 
nc=10 # Number of markets
wpl=2 # Weight of pl in K-means
Random.seed!(12)
R=Clustering.kmeans([zip./1000.0 pl*wpl]',nc);
#R=ParallelKMeans.kmeans([x y pl*wpl]',nc; max_iters=1000);
#R=Clustering.kmeans([x y pl*wpl]',nc);
# R.centers
# R.converged
markets=R.assignments;

##Computing average v and pl per market
vbar=zeros(nc); plbar=zeros(nc);
for i=1:nc
    vbar[i]=mean(v[markets.==i]); plbar[i]=mean(pl[markets.==i]);
end
#Sorting markets based on plbar
marketperm=sortperm(plbar);
markets_sort=zeros(size(markets)); # Markets sorted according to plbar
for m in 1:nc
    markets_sort[markets.==marketperm[m]].=1.0*m
end

for i=1:nc
    vbar[i]=mean(v[markets_sort.==i]); plbar[i]=mean(pl[markets_sort.==i]);
end

#Map of clusters
# scatter(x, y, marker_z=markets_sort, color=:lightrainbow, legend=false)
# scatter(markets_sort, pl, color=:lightrainbow, legend=false)
scatter(markets_sort, v, color=:lightrainbow, legend=false)
#Markets 8 and 9 have very similar upper clusters
#We will demean them, pull them together and estimate the measurement error density


#################################### Deconvolution of mismeasured values #################################### 
### Estimation of the pdf of the measurement error
#Office machine. 58 observations
# findfirst(sort(v[markets_sort.==9]).>39)
# dataforme1=sort(v[markets_sort.==8])[379:end]
# dataforme2=sort(v[markets_sort.==9])[211:end]
# dataforme=sort(vcat(dataforme1 .-mean(dataforme1),dataforme2 .-mean(dataforme2)))
#Home machine
dataforme1=sort(v[markets_sort.==8])[379:end]
dataforme2=sort(v[markets_sort.==9])[211:end]
dataforme=sort(vcat(dataforme1 .-mean(dataforme1),dataforme2 .-mean(dataforme2)))
#####################################
histogram(dataforme,bins=-8:1:8)
kerfunB=Epanechnikov
kerfunK=epanechnikov
#kerfun=Normal
# kerfun=Logistic
# kerfun=SymTriangularDist

dens, sup_length=mixed_density_app(dataforme, kerfunB,kerfunK)
tt=-8:0.01:8
plot(tt,dens.(tt))
### Estimation of v(pl,e)
T=5
npoints=100
τ=0.0001
twostep=true
Theta=Array{Vector{Float64}}(undef,nc)
P=Array{Vector{Float64}}(undef,nc)
model = Model(with_optimizer(KNITRO.Optimizer))
for m in 1:nc
    data=sort(v[markets_sort.==m])
    θ_hat, p_hat=mixture_dist1(data, dens, sup_length, T, npoints, τ)
    if twostep
        θ_hat, p_hat=mixture_dist2(data, dens, vcat(θ_hat, p_hat), model)
    end
    Theta[m]=θ_hat
    P[m]=p_hat
end
#################################### Estimation of proxy function #################################### 
f_po=g_proxy(vbar,plbar,3);
f_po=g_proxy(v,pl,5);
t=minimum(v):maximum(v)
plot(20:60,f_po.(20:60))
##Estimation of output
yo=zeros(T,nc);
po=zeros(nc);
for m in 2:nc
    po[m]=f_po(vbar[m])
    yo[:,m]=Theta[m]./po[m]
end
plot(log.(po),log.(vbar./po))
plot(log.(po[2:end]),log.(yo[:,2:end]'))


m=zeros(nc)
for k in 1:nc
    m[k]=vbar[k]-plbar[k]
end

plot(log.(m[2:end]),log.(yo[:,2:end]'))

