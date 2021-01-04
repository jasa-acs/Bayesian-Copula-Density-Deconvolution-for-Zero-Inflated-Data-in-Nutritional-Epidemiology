
#####################################
### Bayesian Copula Deconvolution ###
#####################################

# Codes accompanying "Bayesian Copula Density Deconvolution for Zero-Inflated Data with Applications in Nutritional Epidemiology" by Sarkar, Pati, Mallick and Carroll.
# Codes written by Abhra Sarkar (abhra.sarkar@utexas.edu), last modified on Dec 15, 2019, in Austin, TX

# The current file is for univariate density deconvolution for variables whose surrogates include continuous positive recalls as well as exact zeros. 
# The method uses normalized mixtures of B-spines to model the density of interest, 
# mixtures of moment-restricted normals to model the density of the measurement errors,
# and mixtures of B-spines to model the conditional variability of the measurement errors.
# The method uses additional latent variables to impute the exact zeros, using mixtures of B-spines to model their means.
# See paper for additional details. 



#############
### Input ###
#############

# While running from within the file 'Bayes_Copula_Decon_MVT.R' that implements the multivariate method, these arguments are read from the original file.
# The method can also be independently implemented for univariate deconvolution using the current file.  

# ws <- surrogates comprising continuous, positive values as well as exact zeros 
# xs.lwr <- lower limit of the range of the variable of interest
# xs.upr <- upper limit of the range of the variable of interest
# mis <- no of surrogates for each subject, must be greater than or equal to 3
# z.xs.max <- number of mixture components allowed in the model for the density of interest
# z.us.max <- number of mixture components allowed in the model for the density of the measurement errors
# K.t <- number of B-spline knots for the variance functions modeling conditional variability of the measurement errors
# simsize <- total num of MCMC iterations
# burnin <- burnin for the MCMC iterations
# show_progress <- if TRUE, shows progress by printing every 100th iteartion number, MUST be set at FALSE while running in parrellel from within 'Bayes_Copula_Decon_MVT.R'
# plot_results <- if TRUE, plots the estimated density of interest, the estimated density of measurement errors, the estimated variance function etc., MUST be set at FALSE while running in parrellel from within 'Bayes_Copula_Decon_MVT.R'



##############
### Output ###
##############

# Output comprises a list of the following variables. 
# While running from within the file 'Bayes_Copula_Decon_MVT.R' that implements the multivariate method, these variables are used as.

# ws <- surrogates with exact zeros imputed by continous values by the model 
# knots <- knot-points for constructing the B-splines bases that model the conditional variability of the measurement errors 
# thetas.xs.density <- estimated coefficients of B-splines bases that model the density of the variable of interest 
# thetas <- estimated coefficients of B-splines bases that model the conditional variability of the measurement errors 
# thetas.eps <- estimated coefficients of B-splines bases that model the mean of the additional latent variables introduced to model exact zeros 
# xs <- estimated subject-specific values of the variable of interest
# xs.eps <- estimated subject-specific means of additional latent surrogates used to model exact zeros
# xs.consumption.days <- estimated subject-specific average values of the variable of interest based on continuous and positive surrogates only
# us <-  estimated subject and replicate-specific values of the measurement errors
# z.us <- mixture component labels for the mixture model for the density of the measurement errors
# pi.us <- mixture component probabilities for the mixture model for the density of the measurement errors
# params.us <- mixture component-specific parameters for the mixture model for the density of the measurement errors


UNIV_DECON_EPISODIC = function(ws, xs.lwr, xs.upr, mis, z.us.max, K.t, simsize, burnin, show_progress=TRUE, plot_results=TRUE)
{
  #################################
  ### Priors and Initial Values ###
  #################################
  
  ### Initialization and prior of xs and us
  n = length(mis)
  N = sum(mis)
  inds = rep(1:n,times=mis)
  inds1 = inds2 = numeric(n)
  inds1[1] = 1
  inds2[1] = inds1[1]+mis[1]-1
  for(ii in 2:n)
  {
    inds1[ii] = inds1[ii-1]+mis[ii-1]
    inds2[ii] = inds1[ii] + mis[ii]-1
  }
  wbars = tapply(ws,inds,"mean")  
  xs = as.vector(wbars)
  xs[xs <= xs.lwr] = xs.lwr+0.1
  xs[xs >= xs.upr] = xs.upr-0.1
  start.xs = xs
  us = ws - rep(xs,times=mis)
  range.start.xs = diff(range(xs))
  s2is = as.vector(tapply(ws,inds,var))
  xs.grid = seq(xs.lwr,xs.upr,length=500)
  xs.grid.length = length(xs.grid)
  d.ordinates.xs = numeric(n)
  
  ### Knot points etc
  P.t = P.mat(K.t+1) 	# penalty matrix
  knots.t = seq(xs.lwr,xs.upr,length=K.t)
  delta.t = knots.t[2]-knots.t[1]
  B.basis.store = B.basis(xs.grid,knots.t)

  ### Prior and initialization of s2t and thetas
  alpha.t = 100
  beta.t = 1
  s2t = 0.01
  optim_results = optim(rep(1,K.t+1), fr, NULL, xs, mis, knots.t, P.t, s2t, us, method = "BFGS")
  thetas = current.thetas = start.thetas = optim_results$par
  thetas[1] = current.thetas[1] = start.thetas[1] = -3
  prop.sig.thetas = make.positive.definite(prop.sig.thetas.fn(thetas,xs,mis,us,s2t,K.t,P.t,knots.t,n))
  
  var.grid = seq(xs.lwr,xs.upr,length=100)
  vars = current.vars = t(B.basis(xs,knots.t)%*%exp(current.thetas))
  B.basis.var.grid.knots.t = B.basis(var.grid,knots.t)

  ### Prior and initialization of thetas.xs.density
  alpha.density.t = 10
  beta.density.t = 1
  s2t.xs.density = 0.1
  kde.density.est = kde(xs,eval.points=xs.grid)$estimate  
  optim_results = optim(rep(1,K.t+1), fr.density, NULL, xs.grid, knots.t, P.t, s2t.xs.density, kde.density.est, method = "BFGS")
  thetas.xs.density = current.thetas.xs.density = optim_results$par
  
  TempMat = abs(matrix(rep(xs,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
  close.ind = apply(TempMat,1,which.min) 
  d.ordinates.xs = current.d.ordinates.xs = B.basis.store[close.ind,]%*%B.basis.density.coeffs(thetas.xs.density,delta.t)
  
  ### Prior and initialization for mixture
  simsize.mh.us = 10
  z.us = 1+numeric(N)
  alpha.us = 0.1
  a.us = 3
  b.us = 1
  sigma0.us = 4
  params.us = matrix(c(0.5,0,1,1),nrow=z.us.max,ncol=4,byrow=T)		# unique values
  pi.us = rep(1/z.us.max,z.us.max)
  d.ordinates.us = matrix(0,nrow=N,ncol=z.us.max)
  
  ### Prior and initialization for latent variables for episodic components
  indices_consumed = which(ws!=0)
  indices_not_consumed = setdiff(1:N,indices_consumed)
  num_proxies_consumed = length(indices_consumed)
  num_proxies_not_consumed = length(indices_not_consumed)
  
  alpha.t.eps = 100
  beta.t.eps = 1
  s2t.eps = 0.01
  ybin = numeric(N)
  ybin[indices_consumed] = 1
  probs.emp = tapply(ybin,inds,sum)/mis
  optim_results = optim(numeric(K.t+1), fr.prob.consumption, NULL, wbars, knots.t, P.t, s2t.eps, probs.emp, method = "BFGS")
  Mu0.thetas.eps = thetas.eps = current.thetas.eps.start = optim_results$par
  Sigma00.thetas.eps = 0.1
  Sigma0.thetas.eps = Sigma00.thetas.eps * diag(K.t+1)
  
  
  ###########################################################
  ### Define Latent and Other Variables for Episodic Data ###
  ###########################################################
  
  xs.eps = newlog(xs)
  ws.eps = rnorm(N,2,1)
  us.eps = rnorm(N)
  
  
  
  #########################
  ### Tuning Parameters ###
  #########################
  
  
  sig.tune.thetas.1 = 0
  sig.tune.thetas.2 = 0.1
  
  
  
  
  
  ###############################
  ### Storage for MCMC Output ###
  ###############################
  
  
  proposed.xs = current.xs = xs
  proposed.xs.eps = current.xs.eps = xs.eps 
  proposed.xs.consumption.days = current.xs.consumption.days = xs.consumption.days = xs
  proposed.us.eps = current.us.eps = us.eps
  proposed.us = current.us = us
  
  xs.star = numeric(N)

  es.grid = seq(-3,3,length=500)
  density.xs.est = numeric(xs.grid.length)
  var.es = numeric(2)
  var.est = numeric(length(var.grid))
  density.es.est = numeric(length(es.grid))
  prob.consumption.est = numeric(xs.grid.length)
  
  
  current.likelihood.eps = proposed.likelihood.eps = matrix(1,2,n)
  temp.proposed.us.likelihood = temp.current.us.likelihood = matrix(1,2,N)
  
  thetas.xs.density.est = thetas.est  = thetas.eps.est = numeric(length(thetas))
  thetas.xs.density.MCMC = matrix(0,nrow=simsize,ncol=length(thetas.xs.density))
  thetas.MCMC = matrix(0,nrow=simsize,ncol=length(thetas))
  thetas.eps.MCMC = matrix(0,nrow=simsize,ncol=length(thetas))
  
  
  
  ##################
  ### Start MCMC ###
  ##################
  
  
  for (iii in 1:simsize)
  {
    if((show_progress=TRUE)&&(iii%%100==0))
    	print(iii)
    
    
    ### Updating thetas.xs.density
    proposed.thetas.xs.density = t(rmvnorm(1,current.thetas.xs.density,(diag(rep(sig.tune.thetas.1,(K.t+1)))+sig.tune.thetas.2*prop.sig.thetas)))
    TempMat = abs(matrix(rep(xs,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
    close.ind = apply(TempMat,1,which.min) 

    current.d.ordinates.xs = B.basis.store[close.ind,]%*%B.basis.density.coeffs(current.thetas.xs.density,delta.t)
    proposed.d.ordinates.xs = B.basis.store[close.ind,]%*%B.basis.density.coeffs(proposed.thetas.xs.density,delta.t)
    
    current.log.prior = - t(current.thetas.xs.density)%*%P.t%*%current.thetas.xs.density/(2*s2t.xs.density)
    proposed.log.prior = - t(proposed.thetas.xs.density)%*%P.t%*%proposed.thetas.xs.density/(2*s2t.xs.density)
    
    current.log.likelihood = sum(log(current.d.ordinates.xs))
    proposed.log.likelihood = sum(log(proposed.d.ordinates.xs))
    
    log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
    if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
    if(log(runif(1))<log.mh.ratio)
    {
      thetas.xs.density = current.thetas.xs.density = proposed.thetas.xs.density
      d.ordinates.xs = current.d.ordinates.xs = proposed.d.ordinates.xs
    }
    
    ### Updating s2t.xs.density
    s2t.xs.density = 1/rgamma(1,shape=alpha.density.t+(K.t+1)/2,rate=beta.density.t+t(thetas.xs.density)%*%P.t%*%thetas.xs.density)
    
    
    ### Updating xs (and us)
    proposed.xs = rtnorm(n,mean=current.xs,sd=0.1,lower=xs.lwr,upper=xs.upr)	

    TempMat = abs(matrix(rep(proposed.xs,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
    close.ind = apply(TempMat,1,which.min) 
    proposed.prior = B.basis.store[close.ind,]%*%B.basis.density.coeffs(thetas.xs.density,delta.t)
    current.prior = d.ordinates.xs
    
    proposed.xs.eps = B.basis.store[close.ind,]%*%thetas.eps
    proposed.xs.consumption.days = proposed.xs/pnorm(proposed.xs.eps)
    proposed.xs.consumption.days[proposed.xs.consumption.days>max(knots.t)] = max(knots.t)-range.start.xs/1000  # For fixing numerical issues
    
    TempMat = abs(matrix(rep(proposed.xs.consumption.days,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
    close.ind = apply(TempMat,1,which.min) 
    proposed.vars = B.basis.store[close.ind,]%*%exp(thetas)

    proposed.us.eps = ws.eps-proposed.xs.eps[inds]		
    temp.current.us.likelihood = dnorm(current.us.eps,mean=0,sd=1)
    temp.proposed.us.likelihood = dnorm(proposed.us.eps,mean=0,sd=1)
    current.likelihood.eps = tapply(temp.current.us.likelihood,inds,"prod")
    proposed.likelihood.eps = tapply(temp.proposed.us.likelihood,inds,"prod")
    
    proposed.us = ws-proposed.xs.consumption.days[inds]
    k.us = max(z.us)	
    temp.current.us.likelihood = fu_mixnorm(current.us,mean=0,sd=sqrt(current.vars[inds]),pi.us[1:k.us],params.us[1:k.us,])
    temp.proposed.us.likelihood = fu_mixnorm(proposed.us,mean=0,sd=sqrt(proposed.vars[inds]),pi.us[1:k.us],params.us[1:k.us,])
    current.likelihood = tapply(temp.current.us.likelihood,inds,"prod")
    proposed.likelihood = tapply(temp.proposed.us.likelihood,inds,"prod")
    
    current.likelihood = as.vector(current.likelihood.eps * current.likelihood)
    proposed.likelihood = as.vector(proposed.likelihood.eps * proposed.likelihood)
    
    mh.ratio = (proposed.prior * proposed.likelihood * dtnorm(current.xs,mean=proposed.xs,sd=0.1,lower=xs.lwr,upper=xs.upr)) / (current.prior * current.likelihood * dtnorm(proposed.xs,mean=current.xs,sd=0.1,lower=xs.lwr,upper=xs.upr))

    mh.ratio[is.nan(mh.ratio)] = 0
    
    u = runif(n)
    inds.to.replace = (1:n)[u<mh.ratio]
    xs[inds.to.replace] = current.xs[inds.to.replace] = proposed.xs[inds.to.replace]
    xs.eps[inds.to.replace] = current.xs.eps[inds.to.replace] = proposed.xs.eps[inds.to.replace]
    xs.consumption.days[inds.to.replace] = current.xs.consumption.days[inds.to.replace] = proposed.xs.consumption.days[inds.to.replace] 
    vars[inds.to.replace] = current.vars[inds.to.replace] = proposed.vars[inds.to.replace]
    
    us.eps = current.us.eps = ws.eps-xs.eps[inds]
    us = current.us = ws - xs.consumption.days[inds]
    
    
    
    ### Updating thetas
    proposed.thetas = t(rmvnorm(1,current.thetas,(diag(rep(sig.tune.thetas.1,(K.t+1)))+sig.tune.thetas.2*prop.sig.thetas)))
    proposed.thetas[1] = 2*proposed.thetas[2]-proposed.thetas[3]
    TempMat = abs(matrix(rep(xs.consumption.days,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
    close.ind = apply(TempMat,1,which.min) 
    proposed.vars = B.basis.store[close.ind,]%*%exp(proposed.thetas)
    
    current.log.prior = - t(current.thetas)%*%P.t%*%current.thetas/(2*s2t)
    proposed.log.prior = - t(proposed.thetas)%*%P.t%*%proposed.thetas/(2*s2t)
    
    temp.current.likelihood = d.restricted.mix.norm(us,mean=rep(0,times=N),sd=sqrt(current.vars[inds]),params.us[z.us,])
    temp.proposed.likelihood = d.restricted.mix.norm(us,mean=rep(0,times=N),sd=sqrt(proposed.vars[inds]),params.us[z.us,])
    
    current.log.likelihood = sum(log(temp.current.likelihood))
    proposed.log.likelihood = sum(log(temp.proposed.likelihood))
    
    log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
    if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
    if(log(runif(1))<log.mh.ratio)
      {
      thetas = current.thetas = proposed.thetas
      vars = current.vars = proposed.vars
      }
    
    ### Updating s2t
    s2t = 1/rgamma(1,shape=alpha.t+(K.t+1)/2,rate=beta.t+t(thetas)%*%P.t%*%thetas)
    
    
    
    if(iii>min(100,burnin/2))
    {
    ### Updating z.us
    for(ii in 1:N)
    {			
      prob.us = pi.us * d.restricted.mix.norm(us[ii],mean=0,sd=sqrt(vars[inds[ii]]),params.us)
      if(sum(prob.us)==0) {prob.us=rep(1/z.us.max,z.us.max)}
      z.us[ii] = sample(1:z.us.max,1,TRUE,prob.us)
    }
    
    ### Updating cluster probabilities
    n.kk.us = tabulate(z.us,nbins=z.us.max)
    pi.us = rdirichlet(1,alpha.us/z.us.max+n.kk.us)	
    
    ### Updating params.us
    k.us = max(z.us)                # Number of clusters
    if(iii>2000) simsize.mh.us = 1
    for(rr in 1:simsize.mh.us)
    {
      for(kk in 1:k.us)
      {
        temp = which(z.us==kk)
        uspool = us[temp]
        varspool = vars[inds[temp]]
        
        proposed.params.us = r.tnorm.proposal.params.restricted.mix.norm(params.us[kk,])
        
        proposed.log.prior = dnorm(proposed.params.us[2],0,sqrt(sigma0.us),log=T) - (a.us+1)*(log(proposed.params.us[3]) + log(proposed.params.us[4])) - b.us*(1/proposed.params.us[3]+1/proposed.params.us[4])
        current.log.prior = dnorm(params.us[kk,2],0,sqrt(sigma0.us),log=T) - (a.us+1)*(log(params.us[kk,3]) + log(params.us[kk,4])) - b.us*(1/params.us[kk,3]+1/params.us[kk,4])
        
        temp.proposed.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),proposed.params.us))
        temp.current.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),params.us[kk,]))
        temp.proposed.log.likelihood[is.infinite(temp.proposed.log.likelihood)] = 0
        temp.current.log.likelihood[is.infinite(temp.current.log.likelihood)] = 0
        proposed.log.likelihood = sum(temp.proposed.log.likelihood)
        current.log.likelihood = sum(temp.current.log.likelihood)
        
        log.acc.prob = proposed.log.prior + proposed.log.likelihood - current.log.prior - current.log.likelihood
        if(log(runif(1))<log.acc.prob)
          params.us[kk,] = proposed.params.us
      }
      if(k.us<z.us.max)
        for(kk in (k.us+1):z.us.max)
          params.us[kk,] = r.proposal.params.restricted.mix.norm(1,1,3,3,3,3,3)	
    }
    }
    
   var.es = var.e.fn(pi.us[1:k.us],params.us[1:k.us,])
   params.us[1:k.us,2] = params.us[1:k.us,2]/sqrt(var.es)
   params.us[1:k.us,3:4] = params.us[1:k.us,3:4]/var.es
   
    
    ### Updating ws.eps 
    xs.star = xs.eps[inds]
    ws.eps[indices_not_consumed] = rtnorm(num_proxies_not_consumed,mean=xs.star[indices_not_consumed],sd=1,lower=-Inf,upper=0)
    ws.eps[indices_consumed] = rtnorm(num_proxies_consumed,mean=xs.star[indices_consumed],sd=1,lower=0,upper=Inf)   
    
    ### Updating thetas.eps
    if(iii>burnin)
      {
      TempMat = abs(matrix(rep(xs,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
      close.ind = apply(TempMat,1,which.min) 
      Mu.thetas.eps = numeric(K.t+1)
      Sigma.thetas.eps = matrix(0,K.t+1,K.t+1)
      for(ii in 1:n)
        {
        for(jj in inds1[ii]:inds2[ii])
          Mu.thetas.eps = Mu.thetas.eps + ws.eps[jj] * B.basis.store[close.ind[ii],]
        Sigma.thetas.eps = Sigma.thetas.eps + mis[ii] * B.basis.store[close.ind[ii],] %*% t(B.basis.store[close.ind[ii],])
        }
      Sigma.thetas.eps = solve(P.t/s2t.eps + solve(Sigma0.thetas.eps) + Sigma.thetas.eps)
      Mu.thetas.eps = Sigma.thetas.eps %*% (solve(Sigma0.thetas.eps)%*%Mu0.thetas.eps + Mu.thetas.eps)
      thetas.eps = t(rmvnorm(1,Mu.thetas.eps,Sigma.thetas.eps))
      
      xs.eps = current.xs.eps = B.basis.store[close.ind,]%*%thetas.eps
      xs.consumption.days = current.xs.consumption.days = xs/pnorm(xs.eps)
      xs.consumption.days[xs.consumption.days>max(knots.t)] = max(knots.t)-range.start.xs/1000  # For fixing numerical issues
  
      TempMat = abs(matrix(rep(xs.consumption.days,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
      close.ind = apply(TempMat,1,which.min) 
      vars = current.vars = B.basis.store[close.ind,]%*%exp(thetas)
      }

    ### Updating s2t.eps
    s2t.eps = 1/rgamma(1,shape=alpha.t.eps+(K.t+1)/2,rate=beta.t.eps+t(thetas.eps)%*%P.t%*%thetas.eps)
   
    ### Updating ws when not observed 
    C2.us = rbinom(N,1,prob=c(params.us[z.us,1],1-params.us[z.us,1]))+1
    mu.us = cbind(params.us[,2]*(1-params.us[,1])/sqrt(params.us[,1]^2+(1-params.us[,1])^2),-params.us[,2]*params.us[,1]/sqrt(params.us[,1]^2+(1-params.us[,1])^2))	
    mu.ws = xs.consumption.days[inds] + sqrt(vars[inds])*mu.us[cbind(z.us,C2.us)]
    sigmasq.ws = vars[inds]*params.us[cbind(z.us,C2.us+2)]
    ws[indices_not_consumed] = rnorm(num_proxies_not_consumed,mean=mu.ws[indices_not_consumed],sd=sqrt(sigmasq.ws[indices_not_consumed]))
    
    
    us.eps = current.us.eps = ws.eps-xs.eps[inds]
    us = current.us = (ws-xs.consumption.days[inds])
    
    
    
    thetas.xs.density.MCMC[iii,] = thetas.xs.density
    thetas.MCMC[iii,] = thetas
    thetas.eps.MCMC[iii,] = thetas.eps
    
    if(iii>burnin)
    {		
      density.xs.est = density.xs.est + B.basis.store%*%B.basis.density.coeffs(thetas.xs.density,delta.t)
      
      prob.consumption.est = prob.consumption.est + pnorm(B.basis.store%*%thetas.eps)

      k.us = max(z.us)
      var.es = var.e.fn(pi.us[1:k.us],params.us[1:k.us,])
      density.es.est = density.es.est + d.scaled.restricted.mix.norm(es.grid,0,1,pi.us[1:k.us],params.us[1:k.us,])
      var.est = var.est + B.basis.var.grid.knots.t %*% exp(thetas) * var.es
      thetas.xs.density.est = thetas.xs.density.est + thetas.xs.density
      thetas.est = thetas.est + log(var.es) + thetas
      thetas.eps.est = thetas.eps.est + thetas.eps
    }
    
  }
  
  
  
  density.xs.est = density.xs.est/(simsize-burnin)
  density.es.est = density.es.est/(simsize-burnin)
  var.est = var.est/(simsize-burnin)
  thetas.xs.density.est = thetas.xs.density.est/(simsize-burnin)
  thetas.est = thetas.est/(simsize-burnin)
  thetas.eps.est = thetas.eps.est/(simsize-burnin)
  prob.consumption.est = prob.consumption.est/(simsize-burnin)
  
  
  if(plot_results==TRUE)
  {
    dev.new()
    par(mfrow=c(2,2))
    plot(xs.grid,density.xs.est,xlab="x",ylab="f(x)",type="l",lty=1,col="green3",lwd=3)
    
    plot(es.grid,density.es.est,xlab="e",ylab="f(e)",type="l",lty=1,col="green3",lwd=3)
    points(es.grid,dnorm(es.grid),type="l",lty=1)
    
    plot(xs,s2is,pch="*",xlab="x",ylab="v(x)")
    points(var.grid,var.est,type="l",lty=1,col="blue",lwd=2)
    points(var.grid,B.basis.var.grid.knots.t%*%exp(thetas.est),type="l",lty=1,col="green3",lwd=2)
    
    plot(xs.grid,1-prob.consumption.est,ylim=c(0,1),type="l",xlab="x",ylab="prob",lwd=3,col="blue")
    points(xs.grid,prob.consumption.est,type="l",lty=2,lwd=3,col="green3")
    title(main="Probability of non-consumption")
    par(mfrow=c(1,1))
  }
  
  return(list(ws=ws,
              knots=knots.t, thetas.xs.density=thetas.xs.density, thetas=thetas.est, thetas.eps=thetas.eps.est,
              xs=xs, xs.eps=xs.eps, xs.consumption.days=xs.consumption.days,
              us=us, z.us=z.us, pi.us=pi.us, params.us=params.us))
}

