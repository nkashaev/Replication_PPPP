## Kernel
epanechnikov(x) =  (abs(x) <= 1.0) ? 0.75 * (1.0-x^2) : 0.0
∇epanechnikov(x) =  (abs(x) <= 1.0) ? -1.5 * x : 0.0
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
    param_p=exp.(vars[pardim+1:end])
    F=sum(log(sum(param_p[j]*(epanechnikov(data[i]-param_pi[j]) - epanechnikov(data[i]-param_pi[end])) for j in 1:pardim-1)+epanechnikov(data[i]-param_pi[end])) for i in 1:n)
    return F
end

function ∇loglike!(G,vars,data)
    pardim=Int((length(vars)+1)/2)
    param_pi=vars[1:pardim]
    param_p=exp.(vars[pardim+1:end])
    for m in 1:pardim-1
        G[m]=-sum(param_p[m]*∇epanechnikov(data[i]-param_pi[m])/sum(param_p[j]*(epanechnikov(data[i]-param_pi[j]) - epanechnikov(data[i]-param_pi[end])) for j in 1:pardim-1)+epanechnikov(data[i]-param_pi[end]) for i in 1:n)
        G[pardim+m]=sum(param_p[m]*(epanechnikov(data[i]-param_pi[m]) - epanechnikov(data[i]-param_pi[end]))/sum(param_p[j]*(epanechnikov(data[i]-param_pi[j]) - epanechnikov(data[i]-param_pi[end])) for j in 1:pardim-1)+epanechnikov(data[i]-param_pi[end]) for i in 1:n)
    end
    G[pardim]=-sum((1.0-sum(param_p))*∇epanechnikov(data[i]-param_pi[pardim])/sum(param_p[j]*(epanechnikov(data[i]-param_pi[j]) - epanechnikov(data[i]-param_pi[end])) for j in 1:pardim-1)+epanechnikov(data[i]-param_pi[end]) for i in 1:n)
end

#Generating matrix A
function MatA(θ,data,dens=epanechnikov)
    n=length(data)
    A=zeros(n,length(θ))
    for i in 1:n, j in 1:length(θ)
        A[i,j]=dens(data[i] - θ[j])
    end    
    return A
end

function opt1(θ,data,verb=0)
    n=length(data)
    npar=length(θ)
    
    model=Model(with_optimizer(SCS.Optimizer))
    set_optimizer_attribute(model, "max_iters", 10000)
    set_optimizer_attribute(model, "verbose", verb)
    @variable(model, logp[1:n]);
    @variable(model, p[1:npar]>=0.0);
    @constraint(model, simplexcon, sum(p)==1.0);

    A=zeros(n,npar)
    for i in 1:n, j in 1:npar
        A[i,j]=epanechnikov(data[i] - θ[j])
    end    
    @constraint(model, con[i = 1:n], [logp[i], 1.0, sum(p[k]*A[i,k] for k=1:npar)] in MOI.ExponentialCone());
    @objective(model, Max,sum(logp))
    optimize!(model)
    return value.(p)
end

function supportadj!(θ_hat,p_hat,data,sup_length=2.0,dens=epanechnikov)
    badobs=data[findall([sum(p_hat[j]*(dens(data[i]-θ_hat[j]) - dens(data[i]-θ_hat[end])) for j in 1:T-1)+dens(data[i]-θ_hat[end]) for i in 1:n].==0)]
    badsupport=~(length(badobs)==0)
    m=0
    while badsupport
        for j in eachindex(badobs)
            Δ=badobs[j] .-θ_hat
            l=argmin(abs.(Δ))
            θ_hat[l]=θ_hat[l]+sign(Δ[l])*1.001*(abs(Δ[l])-sup_length/2.0)
        end
        badobs=data[findall([sum(p_hat[j]*(dens(data[i]-θ_hat[j]) - dens(data[i]-θ_hat[end])) for j in 1:T-1)+dens(data[i]-θ_hat[end]) for i in 1:n].==0)]
        badsupport=~(length(badobs)==0)
        
        m=m+1
        if m>1000
            return println("Too many iterations")
        end
    end
    return θ_hat
end
