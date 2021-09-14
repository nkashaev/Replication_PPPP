tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]*"/Application_Epple/"
resultsdir=rootdir*"Results"
using Pkg
Pkg.activate(".")

using LinearAlgebra, DelimitedFiles
using Distributions, Statistics
using Plots
using Clustering

########################### Functions ##########################################
#include(dirfunc*"/dgp.jl")
EppleData=readdlm(rootdir*"/Pittsburgh_post1995.txt", ',', Float64, '\n',header=true)
Data=EppleData[1]


k=1
markets=zeros(1000)
for i=101:953
    if sum(Data[:,7].==i)>9
        markets[k]=i
        k=k+1
    end
end

markets=markets[markets.>100]
#You can change the location of municipalities here.
#Put .>850 for the right tale; 100.< .<250 for the left; do nothing for all
markets=markets[markets.>850]
#markets=markets[250 .>markets.>100]

indtemp=in.(Data[:,7],[markets])
pl=Data[indtemp,2]
v=Data[indtemp,3]
mun=Data[indtemp,7]

outliers=(pl.<15)
pl=pl[outliers]
mun=mun[outliers]
v=v[outliers]

N=length(pl)
#Computing average v and pl per market
M=length(markets) #Number of markets
vbar=zeros(M)
plbar=zeros(M)
for i=1:M
    vbar[i]=mean(v[mun.==markets[i]])
    plbar[i]=mean(pl[mun.==markets[i]])
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
x=minimum(vbar):(maximum(vbar)-minimum(vbar))/1000:maximum(vbar)
y=zeros(length(x))
for i=1:length(x)
y[i]=price_o(x[i],beta)
end
plot(x,y)
#Inputing implied output level
pobar=zeros(M)
for i=1:M
pobar[i]=price_o(vbar[i],beta)
end
plot(vbar,pobar,seriestype = :scatter)
po=zeros(N)
for i=1:M
po[mun.==markets[i]].=pobar[i]
end

yo=v./po

#Clustering into 2 types
Asign=zeros(N)
for i=1:M
    yot=zeros(1,sum(mun.==markets[i]))
    yot[1,:]=yo[mun.==markets[i]]
    R=kmeans(yot,2)
    Asign[mun.==markets[i]]=R.assignments
end

scatter(mun[mun.>920], yo[mun.>920], marker_z=Asign[mun.>920],
        color=:lightrainbow, legend=false)


# Heterogeneity in output for different markets
plot(mun,yo,seriestype = :scatter)
plot(mun[mun.<200],yo[mun.<200],seriestype = :scatter)
plot(mun[mun.>920],yo[mun.>920],seriestype = :scatter)

plot(mun[mun.==940],yo[mun.==940],seriestype = :scatter)

plot(po[mun.>925],yo[mun.>925],seriestype = :scatter)

# Heterogeneity in price of land for different markets
plot(mun,pl,seriestype = :scatter)
plot(mun[mun.<200],pl[mun.<200],seriestype = :scatter)
plot(mun[mun.>850],pl[mun.>850],seriestype = :scatter)
