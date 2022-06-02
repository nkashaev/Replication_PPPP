## Kernel
epanechnikov(x::Real) =  (abs(x) <= 1.0) ? 0.75 * (1.0 - x^2) : 0.0

## CDFs based on kernel
epanechnikovInvCDF(x) =  (0.0 <= x <= 1.0) ? 2.0*sin(asin(2.0*x - 1)/3.0) : println("out of bounds")

function dgp(pis,p,n,seed)
    Random.seed!(seed)
    return sort(pis[rand(Categorical(p),n)] .+ [epanechnikovInvCDF(rand(n)[i]) for i in 1:n])
end



function mixture_dist2(data, dens, param_ini)
    model = Model(with_optimizer(KNITRO.Optimizer))
    set_silent(model)
    set_optimizer_attribute(model,"xtol",1e-8)
    register(model, :dens, 1, dens; autodiff = true)
    npar=Int(length(param_ini)/2)
    @variable(model, x[1:2*npar] >= 0.0)
    @constraint(model, c1, sum(x[npar+1:end])==1.0)
    set_start_value.(x, param_ini)
    @NLobjective(model, Max, sum(log(sum(x[npar+j]*dens(data[i]-x[j])  for j in 1:npar)) for i in 1:length(data)))
    JuMP.optimize!(model)
    return value.(x)
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
    return mixture_dist2(data, dens, param_ini)
end
