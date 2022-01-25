## Kernel
epanechnikov(x::Real) =  (abs(x) <= 1.0) ? 0.75 * (1.0 - x^2) : 0.0
biweight(x::Real) =  (abs(x) <= 1.0) ? 0.9375 * (1.0 - x^2)^2 : 0.0
triweight(x::Real) =  (abs(x) <= 1.0) ? 1.09375 * (1.0 - x^2)^3 : 0.0

# ∇epanechnikov(x) =  (abs(x) <= 1.0) ? -1.5 * x : 0.0
# npdf(x)=exp(-x^2/2.0)/sqrt(2.0*pi)
## CDFs based on kernel
epanechnikovInvCDF(x) =  (0.0 <= x <= 1.0) ? 2.0*sin(asin(2.0*x - 1)/3.0) : println("out of bounds")

function dgp(pis,p,n,seed)
    Random.seed!(seed)
    return sort(pis[rand(Categorical(p),n)] .+ [epanechnikovInvCDF(rand(n)[i]) for i in 1:n])
end


# #Generating matrix A
# function MatA(θ,data,dens)
#     n=length(data)
#     A=zeros(n,length(θ))
#     for i in 1:n, j in 1:length(θ)
#         A[i,j]=dens(data[i] - θ[j])
#     end    
#     return A
# end



# function supportadj!(θ_hat, dens, sup_length, data, pis_ini)
#     n=length(data);
#     badobs=data[findall([sum(dens(data[i]-θ_hat[j]) for j in 1:length(θ_hat)) for i in 1:n].==0)]
#     badsupport=~(length(badobs)==0)
#     m=0
#     while badsupport
#         for j in eachindex(badobs)
#             Δ=badobs[j] .-θ_hat
#             l=argmin(abs.(Δ))
#             θ_hat[l]=θ_hat[l]+sign(Δ[l])*1.001*(abs(Δ[l])-sup_length/2.0)
#             if length(data[findall([sum(dens(data[i]-θ_hat[j]) for j in 1:length(θ_hat)) for i in 1:n].==0)])==0
#                 break
#             end
#         end
#         badobs=data[findall([sum(dens(data[i]-θ_hat[j]) for j in 1:length(θ_hat)) for i in 1:n].==0)]
#         badsupport=~(length(badobs)==0)        
#         m=m+1
#         if m>10
#             return pis_ini
#         end
#     end
#     return θ_hat
# end

# function mixture_dist1(data, dens, sup_length, T, npoints, τ, pis_ini)
#     ## SQP method
#     θ=[minimum(data)+k*(maximum(data)-minimum(data))/npoints for k in 2:npoints-2]
#     #out = mixSQP(MatA(θ,data,dens), maxiter = 10,verbose = false);
#     out1 = mixSQP(MatA(θ,data,dens),x = ones(length(θ))/length(θ),verbose = false);

#     #Taking solutions with large enough weigths
#     p_raw=out1["x"]
#     θ_raw=θ[p_raw .>τ]
#     p_raw=p_raw[p_raw .>τ]

#     #Clustering θ to T types
#     cl_θ=Clustering.kmeans(θ_raw', minimum([length(θ_raw),T]); weights=p_raw)
#     θ_hat=sort(cl_θ.centers[:]) #Intial guess for θ
#     #Support adjustment
#     θ_hat=supportadj!(θ_hat, dens, sup_length, data, pis_ini)
#     #Recomputing p_hat
#     out2 = mixSQP(MatA(θ_hat,data,dens),x = ones(length(θ_hat))/length(θ_hat),verbose=false); 
#     p_hat=out2["x"] 
    
#     return θ_hat, p_hat
# end

function mixture_dist2(data, dens, param_ini)
    model = Model(with_optimizer(KNITRO.Optimizer))
    # empty!(model)
    set_silent(model)
    set_optimizer_attribute(model,"xtol",1e-6) #Right now it is 1e-4
    # set_optimizer_attribute(model,"print_level",0)
    #model = Model(with_optimizer(KNITRO.Optimizer))
    register(model, :dens, 1, dens; autodiff = true)
    #set_optimizer_attribute(model,"outlev",0)
    npar=Int(length(param_ini)/2)
    @variable(model, x[1:2*npar] >= 0.0)
    @constraint(model, c1, sum(x[npar+1:end])==1.0)
    set_start_value.(x, param_ini)
    @NLobjective(model, Max, sum(log(sum(x[npar+j]*dens(data[i]-x[j])  for j in 1:npar)) for i in 1:length(data)))
    JuMP.optimize!(model)
    return value.(x)[1:npar], value.(x)[npar+1:end]
end

function dens_est(x::Real,kerfK,dataforme,bandw)
    f=0.0
    for i in eachindex(dataforme)
            f=f+ kerfK((x-dataforme[i])/bandw)
    end
    return f/(length(dataforme)*bandw)
end


function mixed_density_app(dataforme, kerfunB, kerfunK)
    h0= rule_of_thumb2(kerfunB,dataforme)
    lscv_res = lscv(kerfunB,dataforme,FFT(); hlb=0.01*h0 )
    bandw = Real(minimizer(lscv_res))
    sup_length=Real(maximum(dataforme)-minimum(dataforme)+2*bandw)
    dens(x::Real)=dens_est(x::Real,kerfunK,dataforme,bandw)
    return dens, sup_length
end

# function clusteringME(data,cluster_distance)
#     clusters=Clustering.dbscan(data', cluster_distance)
#     d=10000.0; 
#     dt=0.0
#     basecluster=0;
#     for i in 1:length(clusters)
#         dt=(maximum(data[clusters[i].core_indices]) - minimum(data[clusters[i].core_indices]))    
#         if dt<d
#             d=dt
#             basecluster=i
#         end 
#     end
    
#     return sort(data[clusters[basecluster].core_indices]-mean(data[clusters[basecluster].core_indices]))
# end

function clusteringME2(data,cluster_distance)
    clusters=Clustering.dbscan(data', cluster_distance)
    nclusters=length(clusters)
    dt=zeros(nclusters)
    for i in 1:nclusters
        dt[i]=(maximum(data[union(clusters[i].boundary_indices,clusters[i].core_indices)]) - minimum(data[union(clusters[i].boundary_indices,clusters[i].core_indices)]))    
    end
    dt[dt.<=0.5].=1000.0 
    d,basecluster=findmin(dt)
    for i in 1:nclusters
        dt[i]=(maximum(data[union(clusters[i].boundary_indices,clusters[i].core_indices)]) - minimum(data[union(clusters[i].boundary_indices,clusters[i].core_indices)]))    
    end
    ntypes=Int.(ceil.(dt/d))
    ntypes[basecluster]=1
    
    pis_ini=zeros(sum(ntypes))
    m=1
    for k in 1:nclusters
        for l in 1:ntypes[k]
            pis_ini[m]=minimum(data[union(clusters[k].boundary_indices,clusters[k].core_indices)])+ dt[k]*(l-0.5)/ntypes[k]    
            m=m+1
        end
    end
    baseind=union(clusters[basecluster].boundary_indices,clusters[basecluster].core_indices)
    return sort(data[baseind] .- mean(data[baseind])), pis_ini
end


# function onesimulation(seed)
#     data=dgp(pis,p,n,seed)
#     dataforme, pis_ini=clusteringME2(data,cluster_distance)
#     dens, sup_length=mixed_density_app(dataforme, kerfunB, kerfunK)
#     if onestep
#         θ_hat1, p_hat1=mixture_dist1(data, dens, sup_length, T, npoints, τ, pis_ini)
#     else
#         θ_hat1=pis
#         p_hat1=p
#     end
    
#     if twostep
#         θ_hat, p_hat=mixture_dist2(data, dens, vcat(θ_hat1, p_hat1))
#     else
#         θ_hat=θ_hat1
#         p_hat=p_hat1
#     end
#     return vcat(θ_hat, p_hat)
# end

function onesimulation(seed)
    data=dgp(pis,p,n,seed)
    dataforme, pis_ini=clusteringME2(data,cluster_distance)
    dens, sup_length=mixed_density_app(dataforme, kerfunB, kerfunK)
    if T==length(pis)
        θ_hat1=pis
        p_hat1=p
    else
        θ_hat1=vcat(pis, [pis[1]+(pis[3]-pis[1])*k/(T-length(pis)+1) for k=1:T-length(pis)])
        p_hat1=vcat(p,zeros(T-length(pis)))
    end    
    param_ini=vcat(θ_hat1, p_hat1)
    θ_hat, p_hat=mixture_dist2(data, dens, param_ini)
    return vcat(θ_hat, p_hat)
end
