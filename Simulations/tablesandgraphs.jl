using Pkg
Pkg.activate(".")


using Distributions, Statistics, Random
using CSV, DataFrames
using Plots
# using LaTeXStrings
# using StatsPlots
# using CategoricalArrays
#################################### Dir ###############################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Simulations/"
dirresults=dirmain*"results"
dirtg=dirmain*"tables and graphs"
#################################### Functions ###############################
function table_simul(kerfunK,N,T,pis)
    Table=zeros(T,2*length(N))
    # dgpN= T==4 ? 2 : 1
    for i in 1:length(N)
        n=N[i]
        output=Matrix(CSV.read(dirresults*"/thetap_DGP_$(n)_$(kerfunK)_$(T).csv", DataFrame))[1:T,:]
        # output=Matrix(CSV.read(dirresults*"/thetap_DGP$(dgpN)_$(n)_$(kerfunK)_$(T).csv", DataFrame))
        println(length(unique(findall(output.<-10)[k][2] for k in 1:length(findall(output.<0.0)))))
        # output=output[:,setdiff(1:size(output,2),baddgp)]
        bpi=mean([sort(output[:,k]).-pis for k in 1:size(output,2)])
        rmsepi=sqrt.(mean([(sort(output[:,k]).-pis).^2 for k in 1:size(output,2)]))
        Table[:,2i-1]=bpi
        Table[:,2i]=rmsepi
    end
    return round.(Table,digits=5)
end
################################### True parameter values
pis=[2.0, 3.0, 3.5, 6.0];
p=[0.2,0.25, 0.25, 0.30];
    
#################################### Correct number of types #################################### 
T=4
N=[500, 1000, 1500,5000]
Table_pi_epanechnikov=table_simul("epanechnikov",N,T,pis)
Table_pi_biweight=table_simul("biweight",N,T,pis)
Table_pi_triweight=table_simul("triweight",N,T,pis)
table_simul("epanechnikov",N,5,pis)

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


