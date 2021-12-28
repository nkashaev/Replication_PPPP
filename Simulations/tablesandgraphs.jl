using Pkg
Pkg.activate(".")


using Distributions, Statistics, Random
using CSV, DataFrames
using Plots
using LaTeXStrings
using StatsPlots
using CategoricalArrays
#################################### Dir ###############################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Simulations/"
dirresults=dirmain*"results"
dirtg=dirmain*"tables and graphs"
################################### True parameter values
pis=[2.0, 3.0, 3.5, 6.0];
p=[0.2,0.25, 0.25, 0.30];
    
#################################### Correct number of types #################################### 
T=4
N=[500, 1000, 1500]
Kern=["epanechnikov","biweight","triweight"]
Table4=zeros(length(Kern),2*length(N))
for i in 1:length(N), j in 1:length(Kern)
    n=N[i]; kerfunK=Kern[j]
    output=Matrix(CSV.read(dirresults*"/thetap_DGP1_$(n)_$(kerfunK)_$(T).csv", DataFrame))
    baddgp=unique(findall(output.==-1000)[k][2] for k in 1:length(findall(output.==-1000)))
    output=output[:,setdiff(1:size(output,2),baddgp)]
    bpi=mean([output[1:length(pis),k].-pis for k in 1:size(output,2)])
    rmsepi=sqrt.(mean([(output[1:length(pis),k].-pis).^2 for k in 1:size(output,2)]))
    Table4[j,2i-1]=bpi[4]
    Table4[j,2i]=rmsepi[4]
end

param=vcat(pis,p)
T=4
N=[500, 1000, 1500,5000]
Kern=["epanechnikov","biweight","triweight"]
Table4=zeros(length(param),2*length(N))
kerfunK="epanechnikov"
for i in 1:length(N)
    n=N[i]
    output=Matrix(CSV.read(dirresults*"/thetap_DGP2_$(n)_$(kerfunK)_$(T).csv", DataFrame))
    baddgp=unique(findall(output.<0.01)[k][2] for k in 1:length(findall(output.<0.01)))
    output=output[:,setdiff(1:size(output,2),baddgp)]
    bpi=mean([output[:,k].-param for k in 1:size(output,2)])
    histogram([output[1,k].-param[1] for k in 1:size(output,2)])
    rmsepi=sqrt.(mean([(output[:,k].-param).^2 for k in 1:size(output,2)]))
    Table4[:,2i-1]=bpi
    Table4[:,2i]=rmsepi
end



#################################### Inccorrect number of types #################################### 
T=4
N=[500, 1000, 1500]
Kern=["epanechnikov","biweight","triweight"]
Table4=zeros(length(Kern),2*length(N))
for i in eachindex(N), j in eachindex(Kern)
    output=Matrix(CSV.read(dirresults*"/thetap_DGP1_$(n)_$(kerfunK)_$(T).csv", DataFrame))
    baddgp=unique(findall(output.==-1000)[k][2] for k in 1:length(findall(output.==-1000)))
    output=output[:,~baddgp]
    Table4[j,2i-1]=mean.([Output[1:length(pis),k].-pis for k in 1:size(output,2)])
    Table4[j,2i]=mean.([(Output[1:length(pis),k].-pis).^2 for k in 1:size(output,2)])
end


#################################### Plotting #################################### 
Table1=zeros(2,5)
Table1[1,:]=hcat(mean(v),median(v),std(v),minimum(v),maximum(v))
Table1[2,:]=hcat(mean(pl),median(pl),std(pl),minimum(pl),maximum(pl))
Table1=DataFrame(Table_sum,[:Mean,:Median,:Std,:Min,:Max])
insertcols!(Table_sum,1,:Variable=>["Value per unit of land","price of land"])
CSV.write(dirtg*"/Table1.csv",  Table1)

Fig6=scatter(pl,v,legend=false, xaxis="Price of land",yaxis="Value per unif of land")
savefig(Fig6,dirtg*"/Fig6.pdf")
Fig7=scatter(market, pl, color=:lightrainbow, legend=false, xaxis="Market",yaxis="Price of land")
savefig(Fig7,dirtg*"/Fig7.pdf")
Fig8=scatter(market, v, color=:lightrainbow, legend=false, xaxis="Market",yaxis="Value per unif of land")
savefig(Fig8,dirtg*"/Fig8.pdf")
