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
dirmain=rootdir*"/Application_Epple/"
dirdata=dirmain*"data"
dirresults=dirmain*"results"
dirtg=dirmain*"tables and graphs"
# #################################### Functions ###############################
# include(dirmain*"application_functions.jl")
# include(dirmain*"mixSQP.jl")

#################################### Data #################################### 
nc=8
data=CSV.read(dirdata*"/data_cleaned_$(nc).csv", DataFrame)
pl=data[:,1]; v=data[:,2]; market=data[:,3]; zip=data[:,4];
output=Matrix(CSV.read(dirresults*"/output_level.csv", DataFrame))
price=Matrix(CSV.read(dirresults*"/output_price.csv", DataFrame))[:]
#################################### Plotting #################################### 
Table1=zeros(2,5)
Table1[1,:]=hcat(mean(v),median(v),std(v),minimum(v),maximum(v))
Table1[2,:]=hcat(mean(pl),median(pl),std(pl),minimum(pl),maximum(pl))
Table_sum=DataFrame(Table1,[:Mean,:Median,:Std,:Min,:Max])
insertcols!(Table_sum,1,:Variable=>["Value per unit of land","price of land"])
CSV.write(dirtg*"/Table1.csv",  Table1)

Fig6=scatter(pl,v,legend=false, xaxis="Price of land",yaxis="Value per unif of land")
savefig(Fig6,dirtg*"/Fig6.pdf")
Fig7=scatter(market, pl, color=:lightrainbow, legend=false, xaxis="Market",yaxis="Price of land")
savefig(Fig7,dirtg*"/Fig7.pdf")
Fig8=scatter(market, v, color=:lightrainbow, legend=false, xaxis="Market",yaxis="Value per unif of land")
savefig(Fig8,dirtg*"/Fig8.pdf")
