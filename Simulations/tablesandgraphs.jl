using Pkg
Pkg.activate(".")


using Distributions, Statistics, Random
using CSV, DataFrames
using Plots
#################################### Dir ###############################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Simulations/"
dirresults=dirmain*"results"
dirtg=dirmain*"tables and graphs"
#################################### Functions ###############################
function table_simul_pi(kerfunK,N,T,pis)
    Table=zeros(T,2*length(N))
    for i in 1:length(N)
        n=N[i]
        output=Matrix(CSV.read(dirresults*"/thetap_DGP_$(n)_$(kerfunK)_$(T).csv", DataFrame))[1:T,:]
        bpi=mean([sort(output[:,k]).-pis for k in 1:size(output,2)])
        rmsepi=sqrt.(mean([(sort(output[:,k]).-pis).^2 for k in 1:size(output,2)]))
        Table[:,i]=bpi
        Table[:,length(N)+i]=rmsepi
    end
    return round.(100.0*Table,digits=2)
end

function table_simul_rho(kerfunK,N,T,p)
    Table=zeros(T-1,2*length(N))
    for i in 1:length(N)
        n=N[i]
        output=Matrix(CSV.read(dirresults*"/thetap_DGP_$(n)_$(kerfunK)_$(T).csv", DataFrame))[T+1:end-1,:]
        bp=mean([output[:,k].-p[1:end-1] for k in 1:size(output,2)])
        rmsep=sqrt.(mean([(output[:,k].-p[1:end-1]).^2 for k in 1:size(output,2)]))
        Table[:,i]=bp
        Table[:,length(N)+i]=rmsep
    end
    return round.(100.0*Table,digits=2)
end
################################### True parameter values
pis=[2.0, 3.0, 3.5, 6.0];
p=[0.2,0.25, 0.25, 0.30];    
#################################### Correct number of types #################################### 
T=4
N=[500, 1000, 1500]
Table_pi_epanechnikov=table_simul_pi("epanechnikov",N,T,pis)
Table_p_epanechnikov=table_simul_rho("epanechnikov",N,T,p)
CSV.write(dirtg*"/Table_pi_epanechnikov.csv", DataFrame(Table_pi_epanechnikov,:auto))
CSV.write(dirtg*"/Table_rho_epanechnikov.csv", DataFrame(Table_p_epanechnikov,:auto))
