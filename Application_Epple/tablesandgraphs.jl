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
nc=10
T=4
data=CSV.read(dirdata*"/data_cleaned_$(nc).csv", DataFrame)
pl=data[:,1]; v=data[:,2]; market=data[:,3]; zip=data[:,4];
output=Matrix(CSV.read(dirresults*"/output_level_$(nc)_$(T).csv", DataFrame))
price=Matrix(CSV.read(dirresults*"/output_price_$(nc)_$(T).csv", DataFrame))[:]
#################################### Plotting #################################### 
Table1=zeros(2,5)
Table1[1,:]=hcat(mean(v),median(v),std(v),minimum(v),maximum(v))
Table1[2,:]=hcat(mean(pl),median(pl),std(pl),minimum(pl),maximum(pl))
Table_sum=DataFrame(Table1,[:Mean,:Median,:Std,:Min,:Max])
insertcols!(Table_sum,1,:Variable=>["Value per unit of land","price of land"])
CSV.write(dirtg*"/Table1_$(nc).csv",  Table_sum)

Fig6=scatter(pl,v,legend=false, xaxis="Price of land",yaxis="Value per unif of land")
savefig(Fig6,dirtg*"/Fig_priceofland_value_$(nc).pdf")
Fig7=scatter(market, pl, color=:lightrainbow, legend=false, xaxis="Market",yaxis="Price of land")
savefig(Fig7,dirtg*"/Fig_priceofland_market_$(nc).pdf")
Fig8=scatter(market, v, color=:lightrainbow, legend=false, xaxis="Market",yaxis="Value per unif of land")
savefig(Fig8,dirtg*"/Fig_value_market_$(nc).pdf")

### Footnote 8
findmax([std(pl[market.==k]) for k in 1:8])

#
# plot(log.(price),log.(output')) 

# [(log.(price[m]),log.(output[4,m]')) for m in 1:nc]