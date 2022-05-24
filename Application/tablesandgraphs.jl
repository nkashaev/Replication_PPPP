using Pkg
Pkg.activate(".")


using Distributions, Statistics, Random
using CSV, DataFrames
using Plots
using LaTeXStrings
using StatsPlots
using CategoricalArrays
using JuMP, KNITRO
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
T=4
data=CSV.read(dirdata*"/data_cleaned_$(nc).csv", DataFrame)
pl=data[:,1]; v=data[:,2]; market=data[:,3]; zip=data[:,4];
output=Matrix(CSV.read(dirresults*"/output_level_$(nc)_$(T).csv", DataFrame))
price=Matrix(CSV.read(dirresults*"/output_price_$(nc)_$(T).csv", DataFrame))[:]
#################################### Plotting #################################### 
Table1=zeros(2,5)
Table1[1,:]=hcat(mean(v),median(v),std(v),minimum(v),maximum(v))
Table1[2,:]=hcat(mean(pl),median(pl),std(pl),minimum(pl),maximum(pl))
Table_sum=DataFrame(round.(Table1,digits=2),[:Mean,:Median,:Std,:Min,:Max])
insertcols!(Table_sum,1,:Variable=>["Value per unit of land","price of land"])
CSV.write(dirtg*"/Table1.csv",  Table_sum)

Fig6=scatter(pl,v,legend=false, xaxis="Price of land",yaxis="Value per unif of land")
savefig(Fig6,dirtg*"/Fig_6.pdf")
Fig7=scatter(market, pl, color=:lightrainbow, legend=false, xaxis="Market",yaxis="Price of land")
savefig(Fig7,dirtg*"/Fig_7.pdf")

function rectangle_from_coords(xb,yb,xt,yt)
    [
        xb  yb
        xt  yb
        xt  yt
        xb  yt
        xb  yb
        NaN NaN
    ]
end

rects=[
    rectangle_from_coords(5.7, 39.5, 6.3, 52)
    rectangle_from_coords(6.7, 41, 7.3, 55)
]

Fig8=scatter(market, v, color=:lightrainbow, legend=false, xaxis="Market",yaxis="Value per unif of land")
plot!(rects[:,1], rects[:,2])
savefig(Fig8,dirtg*"/Fig_8.pdf")

### Footnote 34
(maxstd,maxstdmarket)=round.(findmax([std(pl[market.==k]) for k in 1:8]),digits=2)


#Ploting supply
nmarketstoplot=7 #Change this parameter to use less markets
price_toplot=price[nc-nmarketstoplot+1:nc]
output_toplot=output[:,nc-nmarketstoplot+1:nc]
plot(log.(price_toplot),log.(output_toplot')) 
# The arrays below were used in ploting Figure 9
[(log.(price_toplot[m]),log.(output_toplot[1,m]')) for m in 1:nmarketstoplot]
[(log.(price_toplot[m]),log.(output_toplot[2,m]')) for m in 1:nmarketstoplot]
[(log.(price_toplot[m]),log.(output_toplot[3,m]')) for m in 1:nmarketstoplot]
[(log.(price_toplot[m]),log.(output_toplot[4,m]')) for m in 1:nmarketstoplot]

function average_elasticity(y,x)
    elas=0.0
    for i in 1:length(x)-1
        elas=elas+(y[i+1]-y[i])/(x[i+1]-x[i])
    end
    return elas/(length(x)-1)
end 
[average_elasticity(log.(output_toplot[1,:]'),log.(price_toplot))
average_elasticity(log.(output_toplot[2,:]'),log.(price_toplot))
average_elasticity(log.(output_toplot[3,:]'),log.(price_toplot))
average_elasticity(log.(output_toplot[4,:]'),log.(price_toplot))]

# Monotone supply
model_mon=Model(with_optimizer(KNITRO.Optimizer))
npar=length(output_toplot[1,:]')
@variable(model_mon, x[1:npar] >= 0.0)
@constraint(model_mon, c[i=1:npar-1], x[i]<= x[i+1])    
function monotone_supply(y,model_mon)
    @objective(model_mon, Min, sum((y[i]-x[i])^2 for i in 1:npar))
    JuMP.optimize!(model_mon)
    return value.(x)
end

monoutput=zeros(T,nmarketstoplot)
for i in 1:T
    monoutput[i,:]=monotone_supply(output_toplot[i,:]',model_mon)
end

price_mon=price[nc-nmarketstoplot+1:nc]
# The arrays below were used in ploting Figure 10
[(log.(price_toplot[m]),log.(monoutput[1,m]')) for m in 1:nmarketstoplot]
[(log.(price_toplot[m]),log.(monoutput[2,m]')) for m in 1:nmarketstoplot]
[(log.(price_toplot[m]),log.(monoutput[3,m]')) for m in 1:nmarketstoplot]
[(log.(price_toplot[m]),log.(monoutput[4,m]')) for m in 1:nmarketstoplot]


Table2=round.([average_elasticity(log.(monoutput[1,:]'),log.(price_toplot))
average_elasticity(log.(monoutput[2,:]'),log.(price_toplot))
average_elasticity(log.(monoutput[3,:]'),log.(price_toplot))
average_elasticity(log.(monoutput[4,:]'),log.(price_toplot))], digits=2)
