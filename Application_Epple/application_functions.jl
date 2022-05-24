## Kernel
epanechnikov(x::Real) =  (abs(x) <= 1.0) ? 0.75 * (1.0-x^2) : 0.0
∇epanechnikov(x) =  (abs(x) <= 1.0) ? -1.5 * x : 0.0
npdf(x)=exp(-x^2/2.0)/sqrt(2.0*pi)

# This function (1) uses polynomials to estimate pitilde, 
#               (2) solves the differential equation to obtain the proxy function
function g_proxy(x,y,order=2)
    n=length(y)
    # No sonstant since price of land should be zero if the value is zero 
    #(see also Epple et al. 2010, page 918 and footonote 42)
    X=zeros(n,order)
    for i in 1:size(X,2)
        X[:,i]=x.^i./i
    end
    beta=X\y
    return  t->t^beta[1]*exp(sum((t^(k-1)-1.0)*beta[k]/(k-1) for k in 2:length(beta)))   # Normalization as in Epple et al. 2010 (Proposition 4): y(1)=1 implies g(1)=1 since v=g(v)y(v)  
end


function f_basecluster(data,clusters, extratype=0)
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



#Generating matrix A
function MatA(θ,data,dens)
    n=length(data)
    A=zeros(n,length(θ))
    for i in 1:n, j in 1:length(θ)
        A[i,j]=dens(data[i] - θ[j])
    end    
    return A
end



function supportadj!(θ_hat, dens, sup_length, data)
    n=length(data);
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

function mixture_dist1(data, dens, sup_length, T, npoints, τ)
    ## SQP method
    θ=[minimum(data)+k*(maximum(data)-minimum(data))/npoints for k in 2:npoints-2]
    out = mixSQP(MatA(θ,data,dens), maxiter = 10,verbose = false);
    out1 = mixSQP(MatA(θ,data,dens),x = ones(length(θ))/length(θ),verbose = false);

    #Taking solutions with large enough weigths
    p_raw=out1["x"]
    if sum(p_raw .>=τ)<T
        τ=minimum(sort(p_raw,rev=true)[1:T])
    end
    θ_raw=θ[p_raw .>=τ]
    p_raw=p_raw[p_raw .>=τ]

    #Clustering θ to T types
    cl_θ=Clustering.kmeans(θ_raw', minimum([length(θ_raw),T]); weights=p_raw)
    θ_hat=sort(cl_θ.centers[:]) #Intial guess for θ
    #Support adjustment
    θ_hat=supportadj!(θ_hat, dens, sup_length, data)
    #Recomputing p_hat
    out2 = mixSQP(MatA(θ_hat,data,dens),x = ones(length(θ_hat))/length(θ_hat),verbose=false); 
    p_hat=out2["x"] 
    
    return θ_hat, p_hat
end

function mixture_dist2(data, dens, param_ini, model)
    npar=Int(length(param_ini)/2)
    empty!(model)
    @variable(model, x[1:2*npar] >= 0.0)
    @constraint(model, c1, sum(x[npar+1:end])==1.0)
    set_start_value.(x, param_ini)
    @NLobjective(model, Max, sum(log(sum(x[npar+j]*dens(data[i]-x[j])  for j in 1:npar)) for i in 1:length(data)))
    JuMP.optimize!(model)
    return value.(x)[1:npar], value.(x)[npar+1:end]
end

function mixed_density_app(dataforme, kerfunB, kerfunK)
    lscv_res = lscv(kerfunB,dataforme,FFT())
    bandw = Real(minimizer(lscv_res))
    sup_length=Real(maximum(dataforme)-minimum(dataforme)+2*bandw)
    dens(x::Real)=sum( kerfunK((x+Real(minimum(dataforme)-bandw+sup_length/2-z))/bandw) for z in dataforme)/(length(dataforme)*bandw)
    return dens, sup_length
end


