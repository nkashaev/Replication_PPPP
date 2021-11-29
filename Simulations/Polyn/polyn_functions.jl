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



function loglike(vars,data,dens)
    pardim=Int((length(vars)+1)/2)
    param_pi=vars[1:pardim]
    param_p=exp.(vars[pardim+1:end])
    if sum(param_p)>1.0
        F=-10^10
    else
        F=sum(log(sum(param_p[j]*(dens(data[i]-param_pi[j]) - dens(data[i]-param_pi[end])) for j in 1:pardim-1)+dens(data[i]-param_pi[end])) for i in 1:n)
    end
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
function MatA(θ,data,dens)
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

function supportadj!(θ_hat, dens, sup_length, data)
    badobs=data[findall([sum(dens(data[i]-θ_hat[j]) for j in 1:length(θ_hat)) for i in 1:n].==0)]
    badsupport=~(length(badobs)==0)
    m=0
    while badsupport
        for j in eachindex(badobs)
            Δ=badobs[j] .-θ_hat
            l=argmin(abs.(Δ))
            θ_hat[l]=θ_hat[l]+sign(Δ[l])*1.001*(abs(Δ[l])-sup_length/2.0)
            if length(data[findall([sum(dens(data[i]-θ_hat[j]) for j in 1:length(θ_hat)) for i in 1:n].==0)])==0
                break
            end
        end
        badobs=data[findall([sum(dens(data[i]-θ_hat[j]) for j in 1:length(θ_hat)) for i in 1:n].==0)]
        badsupport=~(length(badobs)==0)        
        m=m+1
        if m>10000
            return println("Too many iterations")
        end
    end
    return θ_hat
end

function mixture_dist(data, dens, sup_length, T, npoints, τ, twostep)
    ## SQP method
    θ=[minimum(data)+k*(maximum(data)-minimum(data))/npoints for k in 2:npoints-2]
    out = mixSQP(MatA(θ,data,dens), maxiter = 10,verbose = false);
    out1 = mixSQP(MatA(θ,data,dens),x = ones(length(θ))/length(θ),verbose = false);

    #Taking solutions with large enough weigths
    p_raw=out1["x"]
    θ_raw=θ[p_raw .>τ]
    p_raw=p_raw[p_raw .>τ]

    #Clustering θ to T types
    cl_θ=kmeans(θ_raw', T; weights=p_raw)
    θ_hat=sort(cl_θ.centers[:]) #Intial guess for θ
    #Support adjustment
    θ_hat=supportadj!(θ_hat, dens, sup_length, data)
    #Recomputing p_hat
    out2 = mixSQP(MatA(θ_hat,data,dens),x = ones(length(θ_hat))/length(θ_hat),verbose=false); 
    p_hat=out2["x"] 
    if twostep
        #Final optimization
        param_ini=vcat(θ_hat,log.(ones(length(θ_hat))[1:end-1]/length(θ_hat))) # Log is taken since we use exp(variable) to guarantee that p>0
        #param_ini=vcat(θ_hat,log.(p_hat)[1:end-1]) # Log is taken since we use exp(variable) to guarantee that p>0
        #opt2 = optimize(vars->-loglike(vars,data), param_ini1, iterations=10^5)
        opt = optimize(vars->-loglike(vars,data,dens), param_ini, iterations=10^5)
        sol=Optim.minimizer(opt)
        #sol2=Optim.minimizer(opt3)

        #Checking if the true parameter value gives a better objective function
        println("The solution is better than the true parameter value --",(-loglike(vcat(pis,log.(p[1:end-1])),data,dens)>Optim.minimum(opt)))
        #Final estimates
        θ_hat=sol[1:length(θ_hat)]
        p_hat=vcat(exp.(sol[length(θ_hat)+1:end]),1.0-sum(exp.(sol[length(θ_hat)+1:end])))
    end
    return θ_hat, p_hat
end

function mixed_density(data, kerfun, cluster_distance)
    clusters=dbscan(data', cluster_distance)
    basecluster, d, ntypes, pis_ini=f_basecluster(clusters,0)
    coredata=data[clusters[basecluster].core_indices];
    lscv_res = lscv(kerfun,coredata,FFT())
    bandw = minimizer(lscv_res)
    fkde = kde(kerfun, bandw, coredata, FFT())
    sup_length=maximum(coredata)-minimum(coredata)+2*bandw
    dens(x)=fkde(x+minimum(coredata)-bandw+sup_length/2)
    return dens, sup_length
end


