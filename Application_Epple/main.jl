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
    t->exp(sum([(t-1.0)^k*beta[k]/k for k in 1:length(beta)]))
    return  t->exp(sum([(t-1.0)^k*beta[k]/k for k in 1:length(beta)]))   # Normalization as in Epple et al. 2010 (Proposition 4): y(1)=1 implies g(1)=1 since v=g(v)y(v)  
end

#################################### Data #################################### 
EppleData=readdlm(dirdata*"/Pittsburgh_post1995.txt", ',', Float64, '\n',header=true);
Data=EppleData[1];
#Removing location outliers 
Noutliersxy= (-80.3.<Data[:,9].<-79.67).*(40.25.<Data[:,10]);
Noutliersplv= (Data[:,2].<1100).*(Data[:,3].<60000);
x=Data[Noutliersxy.*Noutliersplv,9]; y=Data[Noutliersxy.*Noutliersplv,10];
pl=Data[Noutliersxy.*Noutliersplv,2];
v=Data[Noutliersxy.*Noutliersplv,3];
N=length(pl);

scatter(pl,v)

#Preparig markets
R=kmeans([x y pl]',50);
markets=R.assignments;
M=nclusters(R) #Number of markets
#Computing average v and pl per market
vbar=zeros(M); plbar=zeros(M);
for i=1:M
    vbar[i]=mean(v[markets.==i]); plbar[i]=mean(pl[markets.==i]);
end
#Sorting markets based on plbar
marketperm=sortperm(plbar);
markets_sort=zeros(size(markets));
for m in 1:M
    markets_sort[markets.==marketperm[m]].=1.0*m
end

#Map of clusters
scatter(x, y, marker_z=markets_sort, color=:lightrainbow, legend=false)
scatter(markets_sort, pl, color=:lightrainbow, legend=false)
scatter(markets_sort, v, color=:lightrainbow, legend=false)

l1=1
l2=20
scatter(markets_sort[l1 .<=markets_sort .<=l2], v[l1 .<=markets_sort .<=l2], color=:lightrainbow, legend=false)

#################################### Estimation of proxy function #################################### 
f_po=g_proxy(v,pl,5);
plot(sort(v),f_po.(sort(v)))
yo=zeros(size(v));
po=zeros(size(v));
for i in 1:length(yo)
    po[i]=f_po(vbar[markets[i]])
    yo[i]=v[i]/po[i]
end

# #Computing po per market
plot(sort(vbar),f_po.(sort(vbar)))
for i=1:M
    pobar[i]=f_po(vbar[i])
end
l1=1
l2=30
scatter(markets_sort[l1 .<=markets_sort .<=l2], po[l1 .<=markets_sort .<=l2], color=:lightrainbow, legend=false)

# marketperm2=sortperm(po);
# markets_sort2=zeros(size(markets));
# for m in 1:M
#     markets_sort2[markets.==marketperm2[m]].=1.0*m
# end



l1=100
l2=150
scatter(markets_sort[l1 .<=markets_sort .<=l2], yo[l1 .<=markets_sort .<=l2], color=:lightrainbow, legend=false)



#Plot the pticing function
xx=minimum(vbar):(maximum(vbar)-minimum(vbar))/100:maximum(vbar)
yy=zeros(length(xx))
for i=1:length(xx)
    yy[i]=log(f_po(xx[i]))
end
plot(xx,yy)

#Inputing implied output level
yo=v./(f_po.(v))
for i=1:M
pobar[i]=price_o(vbar[i],beta)
end
plot(vbar,pobar,seriestype = :scatter)

po=zeros(N)
for i=1:M
po[markets.==i].=pobar[i]
end

yo=v./po

yobar=zeros(N)
for i=1:M
yobar[markets.==i].=mean(yo[markets.==i])
end

#Clustering into 2 types
Prodt=zeros(N)
for i=1:M
    Pr=kmeans(v[markets.==i]',2)
    Prodt[markets.==i]=Pr.assignments
end
yopr=zeros(N,2)
for i=1:M
    c1=mean(yo[(Prodt.==1).*(markets.==i)])
    c2=mean(yo[(Prodt.==2).*(markets.==i)])
    yopr[markets.==i,:].=hcat(minimum([c1,c2]),maximum([c1,c2]))
end

# Heterogeneity in v for different markets
plot(markets[in.(markets,[1:100])],yo[in.(markets,[1:100])],seriestype = :scatter)
a=0
b=3
plot(po[(a.<po.<b)],v[(a.<po.<b)],seriestype = :scatter)
plot(log.(po[a.<po.<b]),log.(yobar[a.<po.<b]),seriestype = :scatter)
scatter(x, y, marker_z=markets,
        color=:lightrainbow, legend=false)

scatter(po, [yopr[:,1] yopr[:,2]],
                color=:lightrainbow, legend=false)

P=sortslices([po plbar2 yopr],dims=1)
plot(P[:,1],P[:,3:4],seriestype = :scatter)

surface( P[:,1], P[:,2], P[:,3], size=[800,480] )

scatter(P[:,1], P[:,2], P[:,3],
                color=:lightrainbow, legend=false)

histogram(markets)

scatter(markets, yo, marker_z=Prod,
        color=:lightrainbow, legend=false)

scatter(markets[Prod.==2], yo[Prod.==2],
                color=:lightrainbow, legend=false)



                PP=sortslices([yobar markets yo], dims=1)
PPP=sortslices(PP[:,2:3], dims=1)
a=5
plot(PP[PP[:,1].<a,1],(PP[PP[:,1].<a,3]).^1.5, seriestype = :scatter)
