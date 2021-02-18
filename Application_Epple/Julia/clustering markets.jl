
using LinearAlgebra, DelimitedFiles
using Distributions, Statistics
using Plots
using Clustering

########################### Dir ################################################
#office
rootdir="/Users/SSC4044-iMac27/Dropbox/PSU/Research/Production function/Application_Epple/Julia"
#home
rootdir="/Users/nailkashaev/Dropbox/PSU/Research/Production function/Application_Epple/Julia"
#dirfunc=rootdir*"/functions"
dirresults=rootdir*"/results"
########################### Functions ##########################################
#include(dirfunc*"/dgp.jl")
EppleData=readdlm(rootdir*"/Pittsburgh_post1995.txt", ',', Float64, '\n',header=true)
Data=EppleData[1]
#Main Data
coordinates=Data[:,9:10]
x=coordinates[:,1];y=coordinates[:,2];
Noutliersxy= (-80.3.<x.<-79.67).*(40.25.<y)

x=x[Noutliersxy]; y=y[Noutliersxy];
pl=Data[Noutliersxy,2]
v=Data[Noutliersxy,3]
ttr=Data[Noutliersxy,4]
N=length(pl)
#Preparig markets
R=kmeans([x y pl/10]',150)
markets=R.assignments
M=nclusters(R) #Number of markets


scatter(x, y, marker_z=markets, color=:lightrainbow, legend=false)
scatter(markets,v, color=:lightrainbow, legend=false)
scatter(markets, pl, color=:lightrainbow, legend=false)


#Computing average v and pl per market
vbar=zeros(M)
plbar=zeros(M)
plbar2=zeros(N)
for i=1:M
    vbar[i]=mean(v[markets.==i])
    plbar[i]=mean(pl[markets.==i])
    plbar2[markets.==i].=mean(pl[markets.==i])
end

X=hcat(ones(M),vbar,vbar.^2)#,vbar.^3,vbar.^4)#,vbar.^5,vbar.^6)#,vbar.^7,vbar.^8)

beta=inv(X'X)X'plbar

function price_o(v,beta)
#g=[1,v,v^2,v^3,v^4,v^5,v^6]
#f=beta'*g;
#f'/v=beta'*[0,1/v,2,3v^1,4v^2,5v^3,6v^4]
return v^(beta[2])exp(beta[3:end]'*[2v])#,3v^2/2,4v^3/3])#,5v^4/4,6v^5/5,7v^6/6,8v^7/7])
end

#Plot the pticing function
xx=minimum(vbar):(maximum(vbar)-minimum(vbar))/100:maximum(vbar)
yy=zeros(length(xx))
for i=1:length(xx)
yy[i]=price_o(xx[i],beta)
end
plot(xx,yy)
#Inputing implied output level
pobar=zeros(M)
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
