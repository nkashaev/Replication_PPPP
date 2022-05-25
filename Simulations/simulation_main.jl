using Pkg
Pkg.activate(".")
using Distributed
using Statistics
using DataFrames, CSV
#numprocs=4 #Office
numprocs=7 #Home
addprocs(numprocs)
@everywhere using Pkg 
@everywhere Pkg.activate(".")

@everywhere begin 
    using LinearAlgebra, DelimitedFiles
    using Distributions, Statistics, Random
    using Plots
    using Clustering
    using LowRankApprox, Printf
    using JuMP, KNITRO, Ipopt, Optim
    using KDEstimation
    #Dir
    tempdir1=@__DIR__
    rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
    dirmain=rootdir*"/Simulations/"
    #Functions
    include(dirmain*"simulation_functions.jl")
    include(dirmain*"mixSQP.jl")

    # Parameter of DGP
    pis=[2.0, 3.0, 3.5, 6.0];
    p=[0.2,0.25, 0.25, 0.30];    
    cluster_distance=0.3
end
kT=ARGS[1]

if kT=="Epanechnikov"
    @everywhere begin 
        kerfunB=Epanechnikov   
        kerfunK=epanechnikov       
    end
elseif kT=="Biweight"
    @everywhere begin 
        kerfunB=Biweight
        kerfunK=biweight
    end
elseif kT=="Triweight"
    @everywhere begin 
        kerfunB=Triweight
        kerfunK=triweight
    end
else 
    error("Wrong kernel")
end


@everywhere n=$(parse(Int16,ARGS[2]))
@everywhere T=$(parse(Int16,ARGS[3]))    

println("T=",T)
println(kerfunB)
println("n=",n)
println("Ready")

@time Results=pmap(onesimulation,1:1000)
Outcome=zeros(2*T,length(Results))
for j in 1:length(Results)
    Outcome[:,j]=vcat(Results[j], -10000.0*ones(size(Outcome,1)-length(Results[j])))
end

## Saving the output
println("Saving the results")
CSV.write(dirmain*"results/thetap_DGP_$(n)_$(kerfunK)_$(T).csv", DataFrame(Outcome,:auto))
println("Done")


