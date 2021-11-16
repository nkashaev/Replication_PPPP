## Different kernels
# unform(x)=(abs(x) <= 1.0) ? 0.5 : 0.0
# triang(x)=(abs(x) <= 1.0) ? 1 - abs(x) : 0.0
epanechnikov(x) =  (abs(x) <= 1.0) ? 0.75 * (1.0-x^2) : 0.0
# quartic(x) = (abs(x)<=1.0) ? (0.9375 * (1.0 - x^2)^2) : 0.0
# triweight(x) = (abs(x)<=1.0) ? 1.09375 * (1.0-x^2)^3  : 0.0
# tricube(x) = (abs(x)<=1.0) ? 0.8641975308641975 * (1.0-abs(x)^3)^3  : 0.0
gaussian(x) = exp( -0.5 * x^2) * 0.3989422804014327
# cosinus(x) = (abs(x)<=1.0) ?   0.7853981633974483    * cos( 1.5707963267948966 * x ) : 0.0
# logistic(x) = 1/(exp(x)+2.0+exp(-u))
# sigmoid(x) =0.6366197723675814/(exp(x)+exp(-u))

## Mean-zero kernel density estimator
# ydata   -- variable
# sh      -- support parameter/bandwith
# baseker -- kernel function
function m0K_est(ydata, sh, baseker::Function=epanechnikov)
    n=length(ydata)
    ydata=ydata .- sum(ydata)/n # recentering the data around the sample average makes it zero mean
    return sum([baseker((t-ydata[i])*sh) for i in 1:n])sh*/n)
end

function objfun(sh,data)
    f=m0K_est(data[i] - pi[j], sh)
    p[j]*m0K_est(data[i] - pi[j], sh)
end

## Computing rule-of-thumb bandwith for 2 order kernels
# xdata   -- covariates (observation x covariate)
# type    -- rule-of-thumb bandwith: the standard rule-of-thumb, Silverman's rule, Scott's rule
# Requires StatsBase
function fix_bandwith(xdata,type)
h=zeros(size(xdata,2))
mux=mean(xdata,dims=1)
sampsize,dx=size(xdata)
if type=="std rule-of-thumb"
  for k=1:size(xdata,2)
    h[k]=1.06*sampsize^(-1/(dx+4))*min(std(xdata[:,k]),iqr(xdata[:,k])/1.349, mean(abs.(xdata[:,k] .-mux[k]))/1.4826)
  end
elseif type=="silverman"
  for k=1:size(xdata,2)
    h[k]=(4/(dx+2))^(1/(dx+4))*sampsize^(-1/(dx+4))*std(xdata[:,k])
  end
elseif type=="scott"
  for k=1:size(xdata,2)
    h[k]=sampsize^(-1/(dx+4))*std(xdata[:,k])
  end
end
return h
end

## Estimation of conditional probabilities
# choices  -- individual choices (observations x time periods)
# xdata    -- covariates (observation x covariate)
# baseker  -- kernel function
# type     -- type of a bandwith: the standard rule-of-thumb, Silverman's rule, Scott's rule, cross-validation
# sym      -- true if P_ijk is symmetric (i.e. P_ijk=P_ikj=...) (E.g., Frum is stable over time)
# a       -- lowerbound for scale
# b       -- upperbound for scale
# tau     -- stoping criterion: stop if abs(b-a)>tau
# discr   -- set true if the number of types is much smaller than the sample size
# condprob -- 3d array of functions: condprob[i,j,k]=P_ijk(\cdot)
function cond_prob_ker(choices, xdata, baseker, type, sym, a, b, tau, discr)
dP=Int(maximum(choices))
condprob=Array{Function}(undef,dP,dP,dP)
if type!="cross-validation"
    h=fix_bandwith(xdata,type)
    if sym
        for i=1:dP
            for j=i:dP
                for k=j:dP
                    z=unique!(collect(permutations([i,j,k])))
                    ydata=zeros(size(xdata,1))
                    for s in z
                        ydata=ydata+((choices[:,1].==s[1]).*(choices[:,2].==s[2]).*(choices[:,3].==s[3]))/length(z)
                    end
                    lamba=NW_estimator( ydata, xdata, h, baseker)
                    for s in z
                        condprob[s[1],s[2],s[3]]=lamba
                    end
                end
            end
        end
    else
        for i=1:dP, j=1:dP, k=1:dP
            ydata=((choices[:,1].==i).*(choices[:,2].==j).*(choices[:,3].==k))
            condprob[i,j,k]=NW_estimator( ydata, xdata, h, baseker)
        end
    end
else
    if sym
        for i=1:dP
            for j=i:dP
                for k=j:dP
                    z=unique!(collect(permutations([i,j,k])))
                    ydata=zeros(size(xdata,1))
                    for s in z
                        ydata=ydata+((choices[:,1].==s[1]).*(choices[:,2].==s[2]).*(choices[:,3].==s[3]))/length(z)
                    end
                    h=cross_val_h_gr(ydata, xdata, baseker, a, b, tau, discr)
                    lamba=NW_estimator( ydata, xdata, h, baseker)
                    for s in z
                        condprob[s[1],s[2],s[3]]=lamba
                    end
                end
            end
        end
    else
        for i=1:dP, j=1:dP, k=1:dP
            ydata=((choices[:,1].==i).*(choices[:,2].==j).*(choices[:,3].==k))
            h=cross_val_h_gr(ydata, xdata, baseker,type,a , b, tau, discr)
            condprob[i,j,k]=NW_estimator( ydata, xdata, h, baseker)
        end
    end
end
return  condprob
end

## Generating a matrix of CCP for every zip code
# xdata    -- covariates (observation x covariate)
# condprob -- 3d array of functions: condprob[i,j,k]=P_ijk(\cdot)
# sym      -- true if P_ijk is symmetric (i.e. P_ijk=P_ikj=...) (E.g., Frum is stable over time)
function Pijk_data(zipcodes,xdata,condprob,sym)
N,dx=size(xdata)
a=unique(i -> xdata[i,:], 1:N)
dP=size(condprob,1)
# Evaluating estimated conditional probabilities at unique prices
Pijk_p=zeros(dP^3*length(a), 3+1+1)
m=1
if sym
    for i=1:dP
        for j=i:dP
            for k=j:dP
                for n in a
                    t=condprob[i,j,k](xdata[n,:])
                    z=unique!(collect(permutations([i,j,k])))
                    for s in z
                        Pijk_p[m,:]= [s[1] s[2] s[3] t zipcodes[n]]
                        m=m+1
                    end
                end
            end
        end
    end
else
    for i=1:dP, j=1:dP, k=1:dP, n in a
        Pijk_p[m,:]= [i j k condprob[i,j,k](xdata[n,:]) zipcodes[n]]
        m=m+1
    end
end
return Pijk_p
end

##Cross validated bandwith. Optimization is based on the golden ratio rule.
# ydata   -- dependent variable
# xdata   -- covariates (observation x covariate)
# baseker -- kernel function
# a       -- lowerbound for scale
# b       -- upperbound for scale
# tau     -- stoping criterion: stop if abs(b-a)>tau
# discr   -- set true if the number of types is much smaller than the sample size
function cross_val_h_gr(ydata, xdata, baseker, a, b, tau, discr)
phic=(1+sqrt(5))/2
c=b-(b-a)/phic
d=a+(b-a)/phic
#Starting bandwith based on the standard rule-of-thumb bandwith
hb=fix_bandwith(xdata,"std rule-of-thumb")

fa=cros_val_obj(a*hb, ydata, xdata, baseker, discr)
fb=cros_val_obj(b*hb, ydata, xdata, baseker, discr)
fc=cros_val_obj(c*hb, ydata, xdata, baseker, discr)
fd=cros_val_obj(d*hb, ydata, xdata, baseker, discr)
h=hb
indic=0
if fa<fc
    println("decrease a")
    indic=1
    return h,indic
elseif fd>fb
    println("increase b")
    indic=2
    return h,indic
end
while abs(b-a)>tau
    if fc<fd
        b=d
        d=c
        c=b-(b-a)/phic
        fc=cros_val_obj(c*hb, ydata, xdata, baseker, discr)
        fd=cros_val_obj(d*hb, ydata, xdata, baseker, discr)
    else
        a=c
        c=d
        d=a+(b-a)/phic
        fc=cros_val_obj(c*hb, ydata, xdata, baseker, discr)
        fd=cros_val_obj(d*hb, ydata, xdata, baseker, discr)
    end
end
h=hb*(a+b)/2
return h, indic
end

##Cross-validation objective function
# h       -- bandwith
# ydata   -- dependent variable
# xdata   -- covariates (observation x covariate)
# baseker -- kernel function
# discr   -- set true if the number of types is much smaller than the sample size
function cros_val_obj(h, ydata, xdata, baseker, discr)
N,dx=size(xdata)
errsq=0.0
if discr
    v=[vcat(ydata[i], xdata[i,:]) for i=1:N]
    d = Dict{eltype(v), Int}()
    for val in v
        d[val] = get(d, val, 0) + 1
    end
    Ind1=indexin(unique(v),v)
    for k in Ind1
        s=vcat(1:k-1,k+1:N)
        fk=sum([ydata[i]*prod([ baseker((xdata[k,g] - xdata[i,g])/h[g]) for g=1:dx]) for i in s])/sum([prod([ baseker((xdata[k,g] - xdata[i,g])/h[g]) for g=1:dx]) for i in s])
        errsq=errsq+d[vcat(ydata[k], xdata[k,:])]*(ydata[k]-fk)^2/N
    end
else
    for k=1:N
        s=vcat(1:k-1,k+1:N)
        fk=sum([ydata[i]*prod([ baseker((xdata[k,g] - xdata[i,g])/h[g]) for g=1:dx]) for i in s])/sum([prod([ baseker((xdata[k,g] - xdata[i,g])/h[g]) for g=1:dx]) for i in s])
        errsq=errsq+(ydata[k]-fk)^2/N
    end
end
return errsq
end

## Cross validation bandwiths
# choices  -- individual choices (observations x time periods)
# xdata    -- covariates (observation x covariate)
# baseker  -- kernel function
# type     -- type of a bandwith: the standard rule-of-thumb, Silverman's rule, Scott's rule, cross-validation
# sym      -- true if P_ijk is symmetric (i.e. P_ijk=P_ikj=...) (E.g., Frum is stable over time)
# a       -- lowerbound for scale
# b       -- upperbound for scale
# tau     -- stoping criterion: stop if abs(b-a)>tau
# discr   -- set true if the number of types is much smaller than the sample size
function cs_bandwith(choices,xdata,baseker,type, sym, a, b, tau, discr)
dP=Int(maximum(choices))
H=zeros(dP^3,size(xdata,2))
IN=zeros(Int8,dP^3)
m=1
if sym
    for i=1:dP
        for j=i:dP
            for k=j:dP
                z=unique!(collect(permutations([i,j,k])))
                ydata=zeros(size(xdata,1))
                for s in z
                    ydata=ydata+((choices[:,1].==s[1]).*(choices[:,2].==s[2]).*(choices[:,3].==s[3]))/length(z)
                end
                H[m,:],IN[m]=cross_val_h_gr(ydata, xdata, baseker, a, b, tau, discr)
                m=m+1
            end
        end
    end
else
    for i=1:dP, j=1:dP, k=1:dP
        ydata=((choices[:,1].==i).*(choices[:,2].==j).*(choices[:,3].==k))
        H[m,:],IN[m]=cross_val_h_gr(ydata, xdata, baseker, a, b, tau, discr)
        m=m+1
    end
end
return H, IN
end
