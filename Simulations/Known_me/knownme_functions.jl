## Kernel
epanechnikov(x) =  (abs(x) <= 1.0) ? 0.75 * (1.0-x^2) : 0.0
npdf(x)=exp(-x^2/2.0)/sqrt(2.0*pi)
## CDFs based on kernel
epanechnikovInvCDF(x) =  (0.0 <= x <= 1.0) ? 2.0*sin(asin(2.0*x - 1)/3.0) : println("out of bounds")

function dgp(pis,p,n,seed)
    Random.seed!(seed)
    return sort(pis[rand(Categorical(p),n)] .+ [epanechnikovInvCDF(rand(n)[i]) for i in 1:n])
end

function f_basecluster(clusters, extratype=0)
    nclusters=length(clusters)
    dt=zeros(nclusters)
    for i in 1:nclusters
        dt[i]=(maximum(data[clusters[i].core_indices]) - minimum(data[clusters[i].core_indices]))    
    end 
    d,basecluster=findmin(dt)
    ntypes=Int.(ceil.(dt/d)) .+ extratype
    ntypes[basecluster]=1
    
    pis_ini=zeros(sum(ntypes))
    m=1
    for k in 1:nclusters
        for l in 1:ntypes[k]
            pis_ini[m]=minimum(data[clusters[k].core_indices])+ dt[k]*(l-0.5)/ntypes[k]    
            m=m+1
        end
    end
    return basecluster, d, ntypes, pis_ini
end

function loglike(vars,data)
    pardim=Int((length(vars)+1)/2)
    param_pi=vars[1:pardim]
    param_pt=vars[pardim+1:end]
    param_p=vcat(param_pt, 1.0 - sum(param_pt))
    # F=sum([log(sum(param_p[j]*epanechnikov(data[i]-param_pi[j]) for j in 1:pardim)) for i in 1:n])
    F=sum([log(sum(param_p[j]*npdf(data[i]-param_pi[j]) for j in 1:pardim)) for i in 1:n])
    return F-10.0^10*(minimum(param_p)<0.0) 
end
