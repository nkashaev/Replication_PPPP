#Preparing data
using Pkg
Pkg.activate(".")

using LinearAlgebra, DelimitedFiles
using Distributions, Statistics, Random
using Clustering
using CSV, DataFrames
using Plots
#################################### Dir ###############################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("PPPP",tempdir1)[end]]
dirmain=rootdir*"/Application_Epple/"
dirdata=dirmain*"data"
dirresults=dirmain*"results"

#################################### Data #################################### 
EppleData=readdlm(dirdata*"/Pittsburgh_post1995.txt", ',', Float64, '\n',header=true);
Data=EppleData[1];
zipcodedata_raw=readdlm(dirdata*"/US Zip Codes from 2013 Government Data.txt", ',', Float64, '\n',header=true)[1];

#Removing locations, pl, and v outliers. Also removing houses with high pl and value 
Noutliersxy= (-80.3.<Data[:,9].<-79.67).*(40.25.<Data[:,10]);
Noutliersplv= (Data[:,2].<7).*(Data[:,3].<55);
println("The percent of dropped observations is ", 100*(1.0-mean(Noutliersplv))) # percent of dropped observations

x=Data[Noutliersxy.*Noutliersplv,9]; y=Data[Noutliersxy.*Noutliersplv,10]; #coordinates
pl=Data[Noutliersxy.*Noutliersplv,2]; v=Data[Noutliersxy.*Noutliersplv,3]; #price of land and value
N=length(pl);

#Zipcodes
Z=(minimum(y).<=zipcodedata_raw[:,2].<=maximum(y)).*(minimum(x).<=zipcodedata_raw[:,3].<=maximum(x))
zip_coord=zipcodedata_raw[Z,:]
zip_l=zeros(N)
for i in 1:N
    D=(x[i] .- zip_coord[:,3]).^2 .+ (y[i] .- zip_coord[:,2]).^2
    zip_l[i]=zip_coord[argmin(D),1]
end
####### Elbow method to define the number of clusters ######################
K=10
inertiaR=zeros(K)
wpl=10 # Weight of pl in K-means
for k in 1:K
    Random.seed!(1)
    R=Clustering.kmeans([x y pl*wpl]',k);
    for i in 1:length(pl)
        inertiaR[k]=inertiaR[k]+norm([x[i] y[i] pl[i]*wpl]' .- R.centers[:,R.assignments[i]])^2
    end
end
kvar=1:K
FigElbow=plot(kvar,inertiaR[kvar]/10^6,legend=false, xaxis="Number of clusters",yaxis="Within cluster sum of squres x mln")

savefig(FigElbow,dirmain*"/tables and graphs/FigElbow.pdf")
#################################### Preparing Markets #################################### 
nc=8 # Number of markets
wpl=10# Weight of pl in K-means
Random.seed!(1) # Gives the same results in kmeans (12)
R=Clustering.kmeans([x y pl*wpl]',nc);
markets=R.assignments;
##Computing average v and pl per market
vbar=zeros(nc); plbar=zeros(nc);
for i=1:nc
    vbar[i]=mean(v[markets.==i]); plbar[i]=mean(pl[markets.==i]);
end
#Sorting markets based on plbar
marketperm=sortperm(plbar);
markets_sort=zeros(size(markets)); # Markets sorted according to plbar
for m in 1:nc
    markets_sort[markets.==marketperm[m]].=1.0*m
end

Data_cl=hcat(pl,v,markets_sort,zip_l)
#Saving results
CSV.write(dirdata*"/data_cleaned_$(nc).csv", DataFrame(Data_cl,[:pl, :v, :markets, :zip]))
