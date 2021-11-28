using Pkg
Pkg.activate(".")

using LinearAlgebra, DelimitedFiles
using Distributions, Statistics
using Plots
using Clustering

#################################### Dir ###############################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirdata=rootdir*"/Application_Epple/data"
dirresults=rootdir*"/Application_Epple/results"
#################################### Functions ###############################
# This function (1) uses polynomials to estimate pitilde, 
#               (2) solves the differential equation to obtain the proxy function
function g_proxy(x,y,order=2)
    n=length(y)
    # No sonstant since price of land should be zero if the value is zero 
    #(see also Epple et al. 2010, page 918 and footonote 42)
    X=zeros(n,order)
    for i in 1:size(X,2)
        X[:,i]=x.^i
    end
    beta=X\y
    #t->exp(sum([(t-1.0)^k*beta[k]/k for k in 1:length(beta)]))
    return  t->exp(sum([(t-1.0)^k*beta[k]/k for k in 1:length(beta)]))   # Normalization as in Epple et al. 2010 (Proposition 4): y(1)=1 implies g(1)=1 since v=g(v)y(v)  
end

# This function computes the sieve log-likelihood
function sieveL(v,markets_sort,ntypes,sdensity,var)
    nmarkets=unique(markets_sort)
    N=length(v)
    θ=zeros(ntypes,nmarkets)
    π=zeros(ntypes,nmarkets)
    for m in 1:nmarkets
        θ[1:ntypes,m]=var[(m-1)*ntypes+1: m*ntypes]
        t=var[ntypes*nmarkets+(m-1)*(ntypes-1)+1: ntypes*nmarkets+m*(ntypes-1)]
        π[1:ntypes-1,m]=t
        π[ntypes,m]=1.0-sum(t)
    end
    α=var[2*(ntypes-1)*nmarkets+1:end]
    L=0.0
    for i in 1:n
        m=markets_sort[i]
        L=L+sum(π[j,m]*sdensity(v[i]-θ[j,m],α) for j in 1:ntypes)
    end
    return L
end

function sieve_density(base)
    f(x,α)=exp(sum(base[j](x)*α[j] for j in 1:length(base)))s
    g(α)=∫f(x,α)dx
    return x,α->f(x,α)/g(α) 
end
#################################### Data #################################### 
EppleData=readdlm(dirdata*"/Pittsburgh_post1995.txt", ',', Float64, '\n',header=true);
Data=EppleData[1];

#Removing location, pl, and v outliers 
Noutliersxy= (-80.3.<Data[:,9].<-79.67).*(40.25.<Data[:,10]);
Noutliersplv= (Data[:,2].<12).*(Data[:,3].<90);
100*(1.0-mean(Noutliersplv)) # percent of dropped observations

x=Data[Noutliersxy.*Noutliersplv,9]; y=Data[Noutliersxy.*Noutliersplv,10]; #coordinates
pl=Data[Noutliersxy.*Noutliersplv,2]; v=Data[Noutliersxy.*Noutliersplv,3]; #price of land and value
N=length(pl);

scatter(pl,v)

#################################### Preparing Markets #################################### 
nc=30 # Number of markets
wpl=2.0 # Weight of pl in K-means
R=kmeans([x y pl*wpl]',nc);
markets=R.assignments;
M=nclusters(R) #Number of markets
#Computing average v and pl per market
vbar=zeros(M); plbar=zeros(M);
for i=1:M
    vbar[i]=mean(v[markets.==i]); plbar[i]=mean(pl[markets.==i]);
end
#Sorting markets based on plbar
marketperm=sortperm(plbar);
markets_sort=zeros(size(markets)); # Markets sorted according to plbar
for m in 1:M
    markets_sort[markets.==marketperm[m]].=1.0*m
end

#Map of clusters
scatter(x, y, marker_z=markets_sort, color=:lightrainbow, legend=false)
scatter(markets_sort, pl, color=:lightrainbow, legend=false)
scatter(markets_sort, v, color=:lightrainbow, legend=false)

l1=5
l2=20
scatter(markets_sort[l1 .<=markets_sort .<=l2], v[l1 .<=markets_sort .<=l2], color=:lightrainbow, legend=false)
#################################### Deconvolution of mismeasured values #################################### 

#################################### Estimation of proxy function #################################### 
f_po=g_proxy(vbar,plbar,5);
t=minimum(v):maximum(v)
plot(20:60,f_po.(20:60))
yo=zeros(size(v));
po=zeros(size(v));
for i in 1:length(yo)
    po[i]=f_po(vbar[markets[i]])
    yo[i]=v[i]/po[i]
end

l1=10
l2=20
scatter(markets_sort[l1 .<=markets_sort .<=l2], yo[l1 .<=markets_sort .<=l2], color=:lightrainbow, legend=false)

