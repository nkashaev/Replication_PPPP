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
function table_simul(kerfunK,N,T,pis)
    Table=zeros(T,2*length(N))
    for i in 1:length(N)
        n=N[i]
        output=Matrix(CSV.read(dirresults*"/thetap_DGP_$(n)_$(kerfunK)_$(T).csv", DataFrame))[1:T,:]
        # println(length(unique(findall(output.<-10)[k][2] for k in 1:length(findall(output.<0.0)))))
        bpi=mean([sort(output[:,k]).-pis for k in 1:size(output,2)])
        rmsepi=sqrt.(mean([(sort(output[:,k]).-pis).^2 for k in 1:size(output,2)]))
        Table[:,i]=bpi
        Table[:,length(N)+i]=rmsepi
    end
    return round.(100.0*Table,digits=2)
end
################################### True parameter values
pis=[2.0, 3.0, 3.5, 6.0];
p=[0.2,0.25, 0.25, 0.30];    
#################################### Correct number of types #################################### 
T=4
N=[500, 1000, 1500]
Table_pi_epanechnikov=table_simul("epanechnikov",N,T,pis)
CSV.write(dirtg*"/Table_epanechnikov.csv", DataFrame(Table_pi_epanechnikov,:auto))
Table_pi_biweight=table_simul("biweight",N,T,pis)
CSV.write(dirtg*"/Table_biweight.csv", DataFrame(Table_pi_biweight,:auto))
Table_pi_triweight=table_simul("triweight",N,T,pis)
CSV.write(dirtg*"/Table_triweight.csv", DataFrame(Table_pi_triweight,:auto))

