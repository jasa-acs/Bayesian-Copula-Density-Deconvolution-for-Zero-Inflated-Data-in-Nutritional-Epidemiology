  
#####################################
### Bayesian Copula Deconvolution ###
#####################################

# Codes accompanying "Bayesian Copula Density Deconvolution for Zero-Inflated Data with Applications in Nutritional Epidemiology" by Sarkar, Pati, Mallick and Carroll.
# Codes written by Abhra Sarkar (abhra.sarkar@utexas.edu), last modified on Dec 15, 2019, in Austin, TX

# The current file is for multivariate density deconvolution for variables 
# first q (>=0) are episodic with continuous positive surrogates as well as exact zero surrogates
# and the remaining p (>=0) are regular with strictly continuous surrogates.
# The method uses normalized mixtures of B-spines to model the marginal densities of interest of the episodic components, 
# mixtures of truncated-normals with shared atoms to model the marginal densities of interest of the regular components, 
# mixtures of moment-restricted normals to model the density of the measurement errors,
# and mixtures of B-spines to model the conditional variability of the measurement errors.
# The method uses additional latent variables to impute the exact zeros, using mixtures of B-spines to model their means.
# See paper for additional details. 



#############
### Input ###
#############

# filename <- name of a csv file containing the data
# Run_Univariate_Models <- if TRUE runs univariate models from within to generate initial values for the multivariate method, may be set at FALSE for additional runs of the multivariate sampler, default is TRUE
# Plot_Results <- if TRUE, produces plots summarizing the results, default is TRUE
# simsize <- total num of MCMC iterations, default is 5000
# burnin <- burnin for the MCMC iterations, default is 3000
# simsize_univ <- total num of MCMC iterations for the univariate samplers, default is 500, when used independently a default of 5000 is recommended
# burnin_univ <- burnin for the MCMC iterations for the univariate samplers, default is 300, when used independently a default of 3000 is recommended
# K.t <- number of B-spline knots for modeling (a) the marginal densities of interest for components whose surrogates include exact zeros, 
#      (b) the variance functions modeling conditional variability of the measurement errors, and (c) the marginal means of the latent surrogates introduced to model exact zeros, default is 10
# z.x.max <- number of mixture components allowed in the model for the marginal densities of components with continuously measured surrogates, default is 10
# z.u.max <- number of mixture components allowed in the model for the marginal densities of the measurement errors, default is 10
# z.x.max.univ <-  number of mixture components allowed in the model for the marginal densities of components with continuously measured surrogates, default is 5     # MUST be <= z.x.max
# z.u.max.univ <- number of mixture components allowed in the model for the marginal densities of the measurement errors, default is 5                                # MUST be <= z.u.max



##############
### Output ###
##############

# Output is saved in a .RData file and comprises the following objects.

# xs <- estimated subject-specific values of the variable of interest
# knots.t <- knot-points for constructing the B-splines bases used in the model
# x.grid <- grid of values on which the marginal densities of the interest were estimated 
# density.x.est <- estimated marginal densities of the interest evaluated on x.grid
# density.x.est.2D.array <- estimated 2-dimensional joint densities of the interest evaluated on x.grid * x.grid
# e.grid <- grid of values on which the marginal densities of the scaled measurement errors were estimated 
# density.e.est <- estimated marginal densities of the scaled measurement errors evaluated on e.grid
# var.grid <- grid of values on which the variance functions are estimated 
# var.est <- estimated variance functions evaluated on var.grid 
# prob.consumption.est <- estimated probabilities of reporting a positive surrogate evaluated on x.grid
# xs.nrmlzd.est <- estimated subject-specific values of the variable of interest normalized by the last regular component
# x.nrmlzd.grid <- grid of values on which the densities of normalized intakes are estimated
# density.x.nrmlzd.est <- estimated densities of normalized intakes evaluated on x.nrmlzd.grid


###########################
### Additional Comments ###
###########################

# For real data sets, it is recommended that the surrogates be scaled to be unitless and have a common range. 
# This can be done by uncommenting lines 103-105 below.
# Non-linear transformations that do not preserve the additivity of the data generating model are not recommended.
# The common range of values of X is specified by "range.xs". 
# Most real data sets will be highly right skewed. 
# So setting the range at 20 units may bring the effective range to much smaller values. 
# For the EATS dataset analyzed in the main paper, for example, this brings the effective range to about 10 units. 
# For the synthetic datasets analyzed in the main paper, the effective range was about 6 units. 
# The grids on which the densities are to be estimates also have to be set accordingly.
# This is done in lines 123-127 below.
# For the real EATS dataset analyzed in the main paper, we set the lower and upper boundaries of the grid at 0 and 10. 
# For the synthetic datasets analyzed in the main paper, we set the lower and upper boundaries of the grid at 0 and 6.
  
# The current file calls the functions UNIV_DECON_EPISODIC() and UNIV_DECON_REGULAR() that share some arguments with it. 
# The recommended defaults are shown here in the function definition below. 
# To facilitate passing of the shared arguments to the nested functions, they should, however, be set outside in the global environment.
# For illustration, see Bayes_Copula_Decon_MVT_Run.R.


  
Bayes_Copula_Decon_MVT <- function(filename,Run_Univariate_Models=TRUE,Plot_Results=TRUE,Save_Results=TRUE,simsize=5000,burnin=3000,simsize_univ=500,burnin_univ=300,K.t=10,z.x.max=10,z.u.max=10,z.x.max.univ=5,z.u.max.univ=5)
  {

  #################
  ### Load Data ###
  #################
  
  Data <- read.csv(filename,header=TRUE)
  
  
  ###########################
  ### Set Size Parameters ###
  ###########################
  
  var.names <- attributes(Data)$names
  dtot <- length(var.names)
  var.names <- var.names[3:dtot]
  d <- length(var.names)
  ids <- Data$inds
  mis <- as.data.frame(table(ids))$Freq
  n <- length(unique(Data$inds))
  inds <- rep(1:n,times=mis)
  N <- sum(mis)
  
  inds1 = inds2 = numeric(n)
  inds1[1] = 1
  inds2[1] = inds1[1]+mis[1]-1
  for(ii in 2:n)
  {
    inds1[ii] = inds1[ii-1]+mis[ii-1]
    inds2[ii] = inds1[ii] + mis[ii]-1
  }
  
  x.lwr = rep(0,d)
  x.upr = rep(6,d)
  ws <- ws.obs <- t(Data[,3:(2+d)])
  
  #range.xs = 20
  #for(jj in 1:d)
  #{
  #  ws[jj,] <- ws.obs[jj,] <- range.xs*(ws[jj,]-min(ws[jj,]))/(max(ws[jj,])-min(ws[jj,])) 
  #}
  
  # No of episodic components & no of regular components + energy, if included
  row_eps = which(rowSums(ws.obs==0)>0) #rows of episodic recalls (with at least one exact zero)
  q = length(row_eps)     #number of epsodic components
  p = d-q                 #number of regular components
  

  
  
  #################################
  ### Priors and Initial Values ###
  #################################
  
  ### Initialization and prior of xs and us
  inds = rep(1:n,times=mis)
  xs = current.xs = proposed.xs = start.xs = s2is = wbars = matrix(0,d,n)
  xs.consumption.days = current.xs.consumption.days = proposed.xs.consumption.days = xs
  xs.trans = matrix(0,p,n)
  us = matrix(0,d,N)
  x.grid.length = 500
  x.grid = matrix(0,d,x.grid.length)
  for(jj in 1:d)
  {
    xs[jj,] = wbars[jj,] = tapply(ws[jj,],inds,"mean")
    xs[jj,xs[jj,]>x.upr[jj]] = max(xs[jj,xs[jj,]<x.upr[jj]])
    s2is[jj,] = as.vector(tapply(ws[jj,],inds,var))
    x.grid[jj,] = seq(x.lwr[jj],x.upr[jj],length=x.grid.length)
  }
  range.start.xs = apply(t(apply(xs,1,range)),1,diff)
  
  if(p>0)
  {
    ### Prior and initialization for mixture
    alpha.x = 0.1
    
    # Independent Normal Inverse-Gamma
    if(p==1)
    {
      mu0.x = mean(apply(as.matrix(t(xs[(q+1):(q+p),])),1,mean))
      sigmasq0.x = mean(apply(as.matrix(t(xs[(q+1):(q+p),])),1,var))
    }else
    {
      mu0.x = mean(apply(xs[(q+1):(q+p),],1,mean))
      sigmasq0.x = mean(apply(xs[(q+1):(q+p),],1,var))
    }
    a.sigmasq.x = 1
    b.sigmasq.x = 1
    
    z.x.max = 10
    z.x = matrix(0,p,n)
    for(jj in 1:p)
      z.x[jj,] = sample(1:z.x.max,size=n,replace=TRUE,prob=rep(1/z.x.max,z.x.max))
    mu.x = sigmasq.x = matrix(0,p,z.x.max)
    pi.x = matrix(0,p,z.x.max)
    n.kk.x = matrix(0,p,z.x.max)
    prob.x = pi.x
    params.x = array(0,dim=c(p,z.x.max,2))
  }
  
  ### Prior and initialization of s2t and thetas
  alpha.t = 100
  beta.t = 1
  s2t = rep(0.01,d)
  P.t = P.mat(K.t+1) 	# penalty matrix
  knots.t = matrix(0,d,K.t)
  var.es = numeric(d)
  delta.t = numeric(d)
  thetas = current.thetas = proposed.thetas = matrix(0,d,K.t+1)
  vars = current.vars = proposed.vars = matrix(0,d,n)
  prop.sig.thetas = array(0,dim=c(d,K.t+1,K.t+1))
  var.grid = matrix(0,d,100)
  
  indices_consumed = ifelse(ws!=0,1,0)
  indices_not_consumed = ifelse(ws==0,1,0)
  num_proxies_consumed = apply(indices_consumed,1,sum)
  num_proxies_not_consumed = apply(indices_not_consumed,1,sum)
  
  
  ### Prior and initialization for mixture
  z.u.max = 10
  z.u = matrix(1,d,N)
  alpha.u = 0.1
  a.u = 3
  b.u = 1
  sigma0.u = 4
  params.u = array(0,dim=c(d,z.u.max,4))
  for(jj in 1:d)
    params.u[jj,,] = matrix(c(0.5,0,1,1),nrow=z.u.max,ncol=4,byrow=T)		# unique values
  pi.u = matrix(0,d,z.u.max) #matrix(1/z.u.max,d,z.u.max)
  n.kk.u = matrix(0,d,z.u.max)
  prob.u = pi.u
  e.grid = seq(-4,4,length=100)
  density.e.est = matrix(0,d,length(e.grid))
  simsize.mh.u = 10
  
  
  if(q>0)
  {
    ### Prior and initialization for density parameters for episodic components 
    thetas.x.density = current.thetas.x.density = proposed.thetas.x.density = matrix(0,q,K.t+1)
    alpha.density.t = rep(100,q)
    beta.density.t = rep(1,q); beta.density.t[2] = 10
    s2t.x.density = rep(0.1,q)
    
    ### Prior and initialization for latent variables for episodic components
    alpha.t.eps = rep(100,q)
    beta.t.eps = rep(1,q)
    s2t.eps = rep(0.01,q)
    xs.eps = matrix(0,q,n)
    thetas.eps = matrix(0,q,K.t+1)
    Mu0.thetas.eps = matrix(0,q,K.t+1)
    Sigma00.thetas.eps = 0.1
    Sigma0.thetas.eps = Sigma00.thetas.eps * diag(K.t+1)
  }
  
  
  if((q>0)&&(Run_Univariate_Models==TRUE))
  {
    print("Running univariate models for episodic components! Thanks for your patience!")
    clusters <- makeCluster(numCores, type="FORK")
    registerDoParallel(clusters)
    foreach(jj=1:q) %dopar%
    {
      filename = paste("HTR_UNIV_episodic",var.names[jj],"simsize",simsize_univ,"burnin",burnin_univ,"seedno",seedno,sep="_")
      filename = paste(filename,".RData",sep="") 
      univ_results = UNIV_DECON_EPISODIC(ws[jj,],x.lwr[jj],x.upr[jj],mis,z.u.max.univ,K.t,simsize_univ,burnin_univ,FALSE,FALSE) 
      save(file=filename,univ_results)
    }
    stopImplicitCluster()
    print("Done running univariate models for episodic components!")
  }
  for(jj in 1:q)
  {
    if(q==0)
      next
    
    filename = paste("HTR_UNIV_episodic",var.names[jj],"simsize",simsize_univ,"burnin",burnin_univ,"seedno",seedno,sep="_")
    filename = paste(filename,".RData",sep="")
    print(filename)
    load(filename)
    
    ws[jj,] = univ_results$ws
    xs[jj,] = univ_results$xs
    us[jj,] = univ_results$us
    knots.t[jj,] = univ_results$knots
    delta.t[jj] = knots.t[jj,2]-knots.t[jj,1]  
    thetas[jj,] = current.thetas[jj,] = univ_results$thetas
    
    thetas.x.density[jj,] = current.thetas.x.density[jj,] = proposed.thetas.x.density[jj,] = univ_results$thetas.xs.density
    
    z.u[jj,] = univ_results$z.us
    pi.u[jj,1:z.u.max.univ] = univ_results$pi.us
    params.u[jj,1:z.u.max.univ,] = univ_results$params.us
    
    xs.eps[jj,] = univ_results$xs.eps
    xs.consumption.days[jj,] = univ_results$xs.consumption.days
    thetas.eps[jj,] = univ_results$thetas.eps
    Mu0.thetas.eps[jj,] = thetas.eps[jj,]
  }
  
  if((p>0)&&(Run_Univariate_Models==TRUE))
  {
    print("Running univariate models for regular components! Thanks again for your patience!")
    clusters <- makeCluster(numCores, type="FORK")
    registerDoParallel(clusters)
    foreach(jj=(q+1):(q+p))	%dopar%
    {
      filename = paste("HTR_UNIV_regular",var.names[jj],"simsize",simsize_univ,"burnin",burnin_univ,"seedno",seedno,sep="_")
      filename = paste(filename,".RData",sep="") 
      univ_results = UNIV_DECON_REGULAR(ws[jj,],x.lwr[jj],x.upr[jj],mis,z.x.max.univ,z.u.max.univ,K.t,simsize_univ,burnin_univ,FALSE,FALSE) 
      save(file=filename,univ_results)
    }
    stopImplicitCluster()
    print("Done running univariate models for regular components!")
  }
  for(jj in (q+1):(q+p))
  {
    if(p==0)
      next
    
    filename = paste("HTR_UNIV_regular",var.names[jj],"simsize",simsize_univ,"burnin",burnin_univ,"seedno",seedno,sep="_")
    filename = paste(filename,".RData",sep="") 
    print(filename)
    load(filename)
    
    xs[jj,] = univ_results$xs
    us[jj,] = univ_results$us
    knots.t[jj,] = univ_results$knots
    delta.t[jj] = knots.t[jj,2]-knots.t[jj,1]  
    thetas[jj,] = current.thetas[jj,] = univ_results$thetas
    
    z.x[jj-q,] = univ_results$z.xs
    pi.x[jj-q,1:z.x.max.univ] = univ_results$pi.xs
    params.x[jj-q,1:z.x.max.univ,] = univ_results$params.xs
    
    mu.x[jj-q,] = params.x[jj-q,,1]
    sigmasq.x[jj-q,] = params.x[jj-q,,2]
    
    z.u[jj,] = univ_results$z.us
    pi.u[jj,1:z.u.max.univ] = univ_results$pi.us
    params.u[jj,1:z.u.max.univ,] = univ_results$params.us 
    
    xs.consumption.days[jj,] = xs[jj,]
  }
  
  for(jj in 1:d)
  {
    current.xs[jj,] = start.xs[jj,] = xs[jj,]
    xs.consumption.days[jj,xs.consumption.days[jj,]>max(knots.t[jj,])] = max(knots.t[jj,])-range.start.xs[jj]/1000   # For fixing numerical issues
    current.xs.consumption.days[jj,] = xs.consumption.days[jj,]
    temp.vars = B.basis(xs.consumption.days[jj,],knots.t[jj,])%*%exp(thetas[jj,])
    temp.vars[temp.vars==0] = max(temp.vars)
    vars[jj,] = current.vars[jj,] = temp.vars
    prop.sig.thetas[jj,,] = make.positive.definite(prop.sig.thetas.fn(thetas[jj,],xs.consumption.days[jj,],mis,us[jj,],s2t[jj],K.t,P.t,knots.t[jj,],n))
    var.grid[jj,] = seq(min(knots.t[jj,]),max(knots.t[jj,]),length=length(var.grid[jj,]))
  }  
  
  ### Redefine shared component parameters mu.x and sigma.sq.x running a small sampler keeping the xs fixed
  
  if(p>0)
  {
    mu.x.new = numeric(z.x.max)
    mu.x.min = Inf
    mu.x.max = -Inf
    for(jj in 1:p)
    {
      mu.x.min = min(mu.x.min,mu.x[jj,z.x[jj,]])
      mu.x.max = max(mu.x.max,mu.x[jj,z.x[jj,]])
    }
    mu.x.new = seq(mu.x.min,mu.x.max,len=z.x.max)
    for(jj in 1:p)
    {
      TempMat = abs(matrix(mu.x[jj,z.x[jj,]],n,z.x.max)-matrix(rep(mu.x.new,n),n,z.x.max,byrow=T))
      z.x[jj,] = apply(TempMat,1,which.min)
    }
    mu.x = mu.x.new
    sigmasq.x = rep(2*var(as.vector(xs))/z.x.max,z.x.max)
    for(iii in 1:100)
    {
      for(jj in 1:p)
      {
        xs.trans.temp.num = pnorm((xs[q+jj,]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))-pnorm((x.lwr[q+jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))
        xs.trans.temp.den = pnorm((x.upr[q+jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))-pnorm((x.lwr[q+jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))
        xs.trans.temp.den[xs.trans.temp.den<=0] = 0.001  # Numerical Stability
        xs.trans.temp = xs.trans.temp.num/xs.trans.temp.den
        xs.trans.temp[xs.trans.temp>=1] = 0.999  # Numerical Stability
        xs.trans.temp[xs.trans.temp<=0] = 0.001  # Numerical Stability
        xs.trans.temp = mu.x[z.x[jj,]]+sqrt(sigmasq.x[z.x[jj,]])*qnorm(xs.trans.temp)
        xs.trans.temp[xs.trans.temp<(x.lwr[q+jj]-10)] = x.lwr[q+jj] - 10  # Numerical Stability
        xs.trans.temp[xs.trans.temp>(x.upr[q+jj]+10)] = x.upr[q+jj] + 10  # Numerical Stability
        xs.trans[jj,] = xs.trans.temp
      }
      
      for(kk in 1:z.x.max)
      {
        xspool = NULL
        for(jj in 1:p)
        {
          temp = which(z.x[jj,]==kk)
          xspool = c(xspool,xs.trans[jj,temp])
        }
        
        sigmasq.temp = 1/(sum(n.kk.x[,kk])/sigmasq.x[kk] + 1/sigmasq0.x)
        mu.temp = (sum(xspool)/sigmasq.x[kk] + mu0.x/sigmasq0.x) * sigmasq.temp
        mu.x[kk] = rnorm(1,mu.temp,sqrt(sigmasq.temp))	
        
        post.a.sigmasq.x = a.sigmasq.x + length(xspool)/2
        post.b.sigmasq.x = b.sigmasq.x + sum((xspool-mu.x[kk])^2)/2
        sigmasq.x[kk] = 1/rgamma(1,shape=post.a.sigmasq.x,rate=post.b.sigmasq.x)
        
        density.numeric.check = dtnorm(x.grid[jj,],mu.x[kk],sqrt(sigmasq.x[kk]),x.lwr[jj],x.upr[jj])
        if((sum(is.nan(density.numeric.check))>0) || (sum(is.infinite(density.numeric.check))>0))
        {
          mu.x[kk] = mean(wbars)
          sigmasq.x[kk] = 10^10
        }
      }
      
      for(jj in 1:p)
      {
        for(ii in 1:n)
        {
          prob.x = pi.x[jj,] * dtnorm(rep(xs[q+jj,ii],z.x.max),mu.x,sqrt(sigmasq.x),lower=x.lwr[q+jj],upper=x.upr[q+jj])
          prob.x[is.nan(prob.x)]=0; 	prob.x[is.infinite(prob.x)]=max(prob.x[is.finite(prob.x)])   # Numerical Stability
          z.x[jj,ii] = sample(z.x.max,1,TRUE,prob.x)   # New z.x[ii] drawn
        }
        n.kk.x[jj,] = tabulate(z.x[jj,],nbins=z.x.max)
        pi.x[jj,] = rdirichlet(1,alpha.x/z.x.max+n.kk.x[jj,])
      }
    }
  }
  
  
  ### Redefine shared component parameters params.u running a small sampler keeping the us fixed
  
  params.u = matrix(c(0.5,0,1,1),nrow=z.u.max,ncol=4,byrow=T)		# unique values
  for(iii in 1:100)
  {
    for(jj in 1:d)  
    {
      for(rr in 1:simsize.mh.u)
      {
        for(kk in 1:z.u.max)
        {
          uspool = varspool = NULL
          for(jj in 1:d)
          {
            temp = which(z.u[jj,]==kk)
            uspool = c(uspool,us[jj,temp])
            varspool = c(varspool,vars[jj,inds[temp]])
          }
          
          proposed.params.u = r.tnorm.proposal.params.restricted.mix.norm(params.u[kk,])
          
          proposed.log.prior = dnorm(proposed.params.u[2],0,sqrt(sigma0.u),log=T) - (a.u+1)*(log(proposed.params.u[3]) + log(proposed.params.u[4])) - b.u*(1/proposed.params.u[3]+1/proposed.params.u[4])
          current.log.prior = dnorm(params.u[kk,2],0,sqrt(sigma0.u),log=T) - (a.u+1)*(log(params.u[kk,3]) + log(params.u[kk,4])) - b.u*(1/params.u[kk,3]+1/params.u[kk,4])
          
          temp.proposed.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),proposed.params.u))
          temp.current.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),params.u[kk,]))
          temp.proposed.log.likelihood[is.infinite(temp.proposed.log.likelihood)] = 0
          temp.current.log.likelihood[is.infinite(temp.current.log.likelihood)] = 0
          proposed.log.likelihood = sum(temp.proposed.log.likelihood)
          current.log.likelihood = sum(temp.current.log.likelihood)
          
          log.acc.prob = proposed.log.prior + proposed.log.likelihood - current.log.prior - current.log.likelihood
          if(log(runif(1))<log.acc.prob)
            params.u[kk,] = proposed.params.u
        }
      }
      for(ii in 1:N)
      {
        prob.u = pi.u[jj,] * d.restricted.mix.norm(us[jj,ii],mean=0,sd=sqrt(vars[jj,inds[ii]]),params.u)
        prob.u[is.nan(prob.u)] = 0
        if(sum(prob.u)==0)
          prob.u=rep(1/z.u.max,z.u.max)
        z.u[jj,ii] = sample(z.u.max,1,TRUE,prob.u)   # New z.u[ii] drawn
      }
      n.kk.u[jj,] = tabulate(z.u[jj,],nbins=z.u.max)
      pi.u[jj,] = rdirichlet(1,alpha.u/z.u.max+n.kk.u[jj,])
    }
  }
  
  
  
  ### Precompute and store Bspline values on fine grids 
  B.basis.var.grid.knots.t = array(0,dim=c(d,length(var.grid[1,]),K.t+1))
  for(jj in 1:d)
    B.basis.var.grid.knots.t[jj,,] = B.basis(var.grid[jj,],knots.t[jj,])
  B.basis.store = array(0,dim=c(d,x.grid.length,K.t+1))
  for(jj in 1:d)
    B.basis.store[jj,,] = B.basis(x.grid[jj,],knots.t[jj,])
  
  ### Initial density estimates 
  density.x.est = matrix(0,d,length(x.grid[1,]))
  density.e.est = matrix(0,d,length(e.grid))
  if(q>0)
  {
    for(jj in 1:q)
      density.x.est[jj,] = B.basis.store[jj,,]%*%B.basis.density.coeffs(thetas.x.density[jj,],delta.t[jj])
  }
  if(p>0)
  {
    for(jj in 1:p)
    {	
      k.x = max(z.x[jj,])
      for(kk in 1:k.x)
        density.x.est[q+jj,] = density.x.est[q+jj,] + pi.x[jj,kk]*dtnorm(x.grid[q+jj,],mu.x[kk],sqrt(sigmasq.x[kk]),x.lwr[jj],x.upr[jj])
      if(sum(is.nan(density.x.est))>0)
        stop("NAN occurred!")
    }
  }
  
  for(jj in 1:d)
  {
    k.u = max(z.u[jj,])
    density.e.est[jj,] = density.e.est[jj,] + d.scaled.restricted.mix.norm(e.grid,0,1,pi.u[jj,1:k.u],params.u[1:k.u,])
  }
  
  
  
  ### Initialize xs.eps, ws.eps and us.eps
  if(q>0)
  {
    ws.eps = matrix(0,q,N)
    for(jj in 1:q)
    {
      xs.star = xs.eps[jj,inds]
      indices_not_consumed_temp = which(indices_not_consumed[jj,]==1)	
      ws.eps[jj,indices_not_consumed_temp] = rtnorm(num_proxies_not_consumed[jj],mean=xs.star[indices_not_consumed_temp],sd=1,lower=-Inf,upper=0)
      indices_consumed_temp = which(indices_consumed[jj,]==1)	
      ws.eps[jj,indices_consumed_temp] = rtnorm(num_proxies_consumed[jj],mean=xs.star[indices_consumed_temp],sd=1,lower=0,upper=Inf)   
    }
    us.eps = current.us.eps = ws.eps-xs.eps[,inds]
  }
  
  ### Prior and initialization for copula correlation matrices
  yxs = matrix(0,d,n)
  if(q>0)
  {
    for(jj in 1:q)
    {
      temp.xsrs = Fx.norm.B.basis(xs[jj,],knots.t[jj,],thetas.x.density[jj,])
      temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
      temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
      yxs[jj,] = qnorm(temp.xsrs)
    }
  }
  if(p>0)
  {  
    for(jj in 1:p)
    {
      temp.xsrs = Fx_mixtnorm(xs[q+jj,],pi.x[jj,],mu.x,sqrt(sigmasq.x),x.lwr[q+jj],x.upr[q+jj])
      temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
      temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
      yxs[q+jj,] = qnorm(temp.xsrs)
    }
  }
  proposed.yxs = current.yxs = yxs
  yus = matrix(0,d,N)
  for(jj in 1:d)
  {
    temp.usrs = Fu_mixnorm(us[jj,],0,sqrt(vars[jj,inds]),pi.u[jj,],params.u)
    temp.usrs[temp.usrs>=1] = 0.999
    temp.usrs[temp.usrs<=0] = 0.001		
    yus[jj,] = qnorm(temp.usrs)
  }
  proposed.yus = current.yus = yus
  
  bx = current.bx = rep(0,d-1)
  thetax = current.thetax = rep(0,(d-1)*(d-2)/2)
  bxset = seq(-0.99,0.99,len=41)
  thetaxset = seq(-3.14,3.14,len=41)
  
  Corr.xs = current.Corr.xs = formCorr.xs(current.bx,current.thetax)
  current.inv.Corr.xs = solve(current.Corr.xs)
  
  be = current.be = rep(0,d-1)
  thetae = current.thetae = rep(0,(d-1)*(d-2)/2)
  beset = seq(-0.99,0.99,len=41)
  thetaeset = seq(-3.14,3.14,len=41)
  
  Corr.es = current.Corr.es = formCorr.es(current.be,current.thetae)
  current.inv.Corr.es = solve(current.Corr.es)
  
  
  
  ###############################
  ### Storage for MCMC Output ###
  ###############################
  
  sig.tune.thetas.1 = 0
  sig.tune.thetas.2 = 1
  
  density.x.est = matrix(0,d,length(x.grid[1,]))
  density.e.est = matrix(0,d,length(e.grid))
  var.true = var.est = matrix(0,d,length(var.grid[1,]))
  var.e = numeric(d)
  Corr.xs.est = matrix(0,d,d)
  Corr.es.est = matrix(0,d,d)
  
  proposed.xs = current.xs = xs
  proposed.us = current.us = us
  
  proposed.prior = current.prior = matrix(1,d,n)
  current.copula.prior = proposed.copula.prior = rep(1,n)
  proposal.p.to.c = proposal.c.to.p = matrix(1,d,n)
  current.marginal.likelihood = proposed.marginal.likelihood = matrix(1,d,n)
  current.copula.likelihood = proposed.copula.likelihood = matrix(1,d,n)
  d.ordinates.x = current.d.ordinates.x = proposed.d.ordinates.x = matrix(0,q,n)
  if(q>0)
  {
    proposed.xs.eps = current.xs.eps = xs.eps
    proposed.us.eps = current.us.eps = us.eps
    current.marginal.likelihood.eps = proposed.marginal.likelihood.eps = matrix(1,d,n)
    prob.consumption.est = matrix(0,q,x.grid.length)
    for(jj in 1:q)
    {
      TempMat = abs(matrix(rep(xs[jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[jj,],n),n,x.grid.length,byrow=T))
      close.ind = apply(TempMat,1,which.min) 
      d.ordinates.x[jj,] = B.basis.store[jj,close.ind,]%*%B.basis.density.coeffs(thetas.x.density[jj,],delta.t[jj])
    }
  }
  
  
  ##################
  ### Start MCMC ###
  ##################
  
  
  
  for (iii in 1:simsize)
  {
    if(iii%%100==0)
      print(iii)
    
    if(q>0)
    {      
      for(jj in 1:q)
      {
        ### Updating thetas.xs.density
        
        proposed.thetas.x.density[jj,] = t(rmvnorm(1,current.thetas.x.density[jj,],(diag(rep(sig.tune.thetas.1,(K.t+1)))+sig.tune.thetas.2*prop.sig.thetas[jj,,])))
        TempMat = abs(matrix(rep(xs[jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[jj,],n),n,x.grid.length,byrow=T))
        close.ind = apply(TempMat,1,which.min) 
        
        current.d.ordinates.x[jj,] = B.basis.store[jj,close.ind,]%*%B.basis.density.coeffs(current.thetas.x.density[jj,],delta.t[jj])
        proposed.d.ordinates.x[jj,] = B.basis.store[jj,close.ind,]%*%B.basis.density.coeffs(proposed.thetas.x.density[jj,],delta.t[jj])
        
        current.log.prior = - t(current.thetas.x.density[jj,])%*%P.t%*%current.thetas.x.density[jj,]/(2*s2t.x.density[jj])
        proposed.log.prior = - t(proposed.thetas.x.density[jj,])%*%P.t%*%proposed.thetas.x.density[jj,]/(2*s2t.x.density[jj])
        
        current.log.likelihood = sum(log(current.d.ordinates.x[jj,]))
        proposed.log.likelihood = sum(log(proposed.d.ordinates.x[jj,]))
        
        log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
        if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
        if(log(runif(1))<log.mh.ratio)
        {
          thetas.x.density[jj,] = current.thetas.x.density[jj,] = proposed.thetas.x.density[jj,]
          d.ordinates.x[jj,] = current.d.ordinates.x[jj,] = proposed.d.ordinates.x[jj,]
          
          ### Updating yxs
          
          temp.xsrs = Fx.norm.B.basis(xs[jj,],knots.t[jj,],thetas.x.density[jj,])
          temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
          temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
          yxs[jj,] = current.yxs[jj,] = qnorm(temp.xsrs)
        }      
        
        
        ### Updating s2t.x.density
        
        s2t.x.density[jj] = 1/rgamma(1,shape=alpha.density.t[jj]+(K.t+1)/2,rate=beta.density.t[jj]+t(thetas.x.density[jj,])%*%P.t%*%thetas.x.density[jj,])
      }
    }
    
    
    if(p>0)
    {      
      ### Updating mu.x, sigmasq.x
      
      for(kk in 1:z.x.max)
      {
        xspool = NULL
        for(jj in 1:p)
        {
          temp = which(z.x[jj,]==kk)
          xspool = c(xspool,xs.trans[jj,temp])
        }
        
        sigmasq.temp = 1/(sum(n.kk.x[,kk])/sigmasq.x[kk] + 1/sigmasq0.x)
        mu.temp = (sum(xspool)/sigmasq.x[kk] + mu0.x/sigmasq0.x) * sigmasq.temp
        mu.x[kk] = rnorm(1,mu.temp,sqrt(sigmasq.temp))	
        
        post.a.sigmasq.x = a.sigmasq.x + length(xspool)/2
        post.b.sigmasq.x = b.sigmasq.x + sum((xspool-mu.x[kk])^2)/2
        sigmasq.x[kk] = 1/rgamma(1,shape=post.a.sigmasq.x,rate=post.b.sigmasq.x)
        
        density.numeric.check = dtnorm(x.grid[jj,],mu.x[kk],sqrt(sigmasq.x[kk]),x.lwr[jj],x.upr[jj])
        if((sum(is.nan(density.numeric.check))>0) || (sum(is.infinite(density.numeric.check))>0))
        {
          mu.x[kk] = mean(wbars)
          sigmasq.x[kk] = 10^10
        }
      }
      
      for(jj in 1:p)
      {
        ### Updating z.x and pi.x
        
        for(ii in 1:n)
        {
          prob.x = pi.x[jj,] * dtnorm(rep(xs[q+jj,ii],z.x.max),mu.x,sqrt(sigmasq.x),lower=x.lwr[q+jj],upper=x.upr[q+jj])
          prob.x[is.nan(prob.x)]=0; 	prob.x[is.infinite(prob.x)]=max(prob.x[is.finite(prob.x)])   # Numerical Stability
          z.x[jj,ii] = sample(z.x.max,1,TRUE,prob.x)
        }
        n.kk.x[jj,] = tabulate(z.x[jj,],nbins=z.x.max)
        pi.x[jj,] = rdirichlet(1,alpha.x/z.x.max+n.kk.x[jj,])	
        
        
        ### Updating yxs
        
        temp.xsrs = Fx_mixtnorm(xs[q+jj,],pi.x[jj,],mu.x,sqrt(sigmasq.x),x.lwr[q+jj],x.upr[q+jj])
        temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
        temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
        yxs[q+jj,] = current.yxs[q+jj,] = qnorm(temp.xsrs)
      }
    }
    
    
    
    ### Updating params.u
    
    k.u = max(z.u)                # Number of clusters
    if(iii>2000)
      simsize.mh.u = 1
    for(rr in 1:simsize.mh.u)
    {
      for(kk in 1:k.u)
      {
        uspool = varspool = NULL
        for(jj in 1:d)
        {
          temp = which(z.u[jj,]==kk)
          uspool = c(uspool,us[jj,temp])
          varspool = c(varspool,vars[jj,inds[temp]])
        }
        
        proposed.params.u = r.tnorm.proposal.params.restricted.mix.norm(params.u[kk,])
        
        proposed.log.prior = dnorm(proposed.params.u[2],0,sqrt(sigma0.u),log=T) - (a.u+1)*(log(proposed.params.u[3]) + log(proposed.params.u[4])) - b.u*(1/proposed.params.u[3]+1/proposed.params.u[4])
        current.log.prior = dnorm(params.u[kk,2],0,sqrt(sigma0.u),log=T) - (a.u+1)*(log(params.u[kk,3]) + log(params.u[kk,4])) - b.u*(1/params.u[kk,3]+1/params.u[kk,4])
        
        temp.proposed.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),proposed.params.u))
        temp.current.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),params.u[kk,]))
        temp.proposed.log.likelihood[is.infinite(temp.proposed.log.likelihood)] = 0
        temp.current.log.likelihood[is.infinite(temp.current.log.likelihood)] = 0
        proposed.log.likelihood = sum(temp.proposed.log.likelihood)
        current.log.likelihood = sum(temp.current.log.likelihood)
        
        log.acc.prob = proposed.log.prior + proposed.log.likelihood - current.log.prior - current.log.likelihood
        if(log(runif(1))<log.acc.prob)
          params.u[kk,] = proposed.params.u
      }
    }
    if(k.u<z.u.max)
      for(kk in (k.u+1):z.u.max)
        params.u[kk,] = r.proposal.params.restricted.mix.norm(1,1,sigma0.u,a.u,b.u,a.u,b.u)	
    
    
    for(jj in 1:d)
    {
      ### Updating z.u
      
      for(ii in 1:N)
      {
        prob.u = pi.u[jj,] * d.restricted.mix.norm(us[jj,ii],mean=0,sd=sqrt(vars[jj,inds[ii]]),params.u)
        prob.u[is.nan(prob.u)] = 0
        if(sum(prob.u)==0)
          prob.u=rep(1/z.u.max,z.u.max)
        z.u[jj,ii] = sample(z.u.max,1,TRUE,prob.u)
      }
      n.kk.u[jj,] = tabulate(z.u[jj,],nbins=z.u.max)
      pi.u[jj,] = rdirichlet(1,alpha.u/z.u.max+n.kk.u[jj,])
      
      
      ### Updating yus
      
      temp.usrs = Fu_mixnorm(us[jj,],0,sqrt(vars[jj,inds]),pi.u[jj,],params.u)
      temp.usrs[temp.usrs>=1] = 0.999
      temp.usrs[temp.usrs<=0] = 0.001		
      yus[jj,] = current.yus[jj,] = qnorm(temp.usrs)
    }
    
    
    
    ### Updating the Gaussian copula correlation parameters
    
    for(tt in 1:(d-1))
    {
      proposed.bx = current.bx
      proposed.bx[tt] = sample(c(current.bx[tt],current.bx[tt]+c(-1,1)*2*0.99/40),size=1)
      while(abs(proposed.bx[tt])>0.99)
        proposed.bx[tt] = sample(c(current.bx[tt],current.bx[tt]+c(-1,1)*2*0.99/40),size=1)		
      proposed.Corr.xs = formCorr.xs(proposed.bx,thetax)
      proposed.inv.Corr.xs = solve(proposed.Corr.xs)
      
      proposed.log.likelihood = -(n/2)*log(1-proposed.bx[tt]*proposed.bx[tt])
      current.log.likelihood = -(n/2)*log(1-current.bx[tt]*current.bx[tt])
      for(ii in 1:n)
      {
        current.log.likelihood = current.log.likelihood - 0.5 * t(yxs[,ii]) %*% current.inv.Corr.xs %*% yxs[,ii]
        proposed.log.likelihood = proposed.log.likelihood - 0.5 * t(yxs[,ii]) %*% proposed.inv.Corr.xs %*% yxs[,ii]
      }
      log.mh.ratio = proposed.log.likelihood - current.log.likelihood	
      if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
      if(log(runif(1))<log.mh.ratio)
      {
        bx[tt] = current.bx[tt] = proposed.bx[tt]
        Corr.xs = current.Corr.xs = proposed.Corr.xs
        inv.Corr.xs = current.inv.Corr.xs = proposed.inv.Corr.xs
      }	
    }
    
    for(tt in 1:((d-2)*(d-1)/2))
    {
      proposed.thetax = current.thetax
      proposed.thetax[tt] = sample(c(current.thetax[tt],current.thetax[tt]+c(-1,1)*2*3.14/40),size=1)
      while(abs(proposed.thetax[tt])>3.14)
        proposed.thetax[tt] = sample(c(current.thetax[tt],current.thetax[tt]+c(-1,1)*2*0.99/40),size=1)		
      proposed.Corr.xs = formCorr.xs(bx,proposed.thetax)
      proposed.inv.Corr.xs = solve(proposed.Corr.xs)
      
      proposed.log.likelihood = 0		
      current.log.likelihood = 0
      for(ii in 1:n)
      {
        current.log.likelihood = current.log.likelihood - 0.5 * t(yxs[,ii]) %*% current.inv.Corr.xs %*% yxs[,ii]
        proposed.log.likelihood = proposed.log.likelihood - 0.5 * t(yxs[,ii]) %*% proposed.inv.Corr.xs %*% yxs[,ii]
      }
      log.mh.ratio = proposed.log.likelihood - current.log.likelihood	
      if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
      if(log(runif(1))<log.mh.ratio)
      {
        thetax[tt] = current.thetax[tt] = proposed.thetax[tt]
        Corr.xs = current.Corr.xs = proposed.Corr.xs
        inv.Corr.xs = current.inv.Corr.xs = proposed.inv.Corr.xs
      }	
    }
    
    
    for(tt in 1:(d-1))
    {
      proposed.be = current.be
      proposed.be[tt] = sample(c(current.be[tt],current.be[tt]+c(-1,1)*2*0.99/40),size=1)
      while(abs(proposed.be[tt])>0.99)
        proposed.be[tt] = sample(c(current.be[tt],current.be[tt]+c(-1,1)*2*0.99/40),size=1)		
      proposed.Corr.es = formCorr.es(proposed.be,thetae)
      proposed.inv.Corr.es = solve(proposed.Corr.es)
      
      proposed.log.likelihood = -(N/2)*log(1-proposed.be[tt]*proposed.be[tt])		
      current.log.likelihood = -(N/2)*log(1-current.be[tt]*current.be[tt])
      for(ii in 1:N)
      {
        current.log.likelihood = current.log.likelihood - 0.5 * t(yus[,ii]) %*% current.inv.Corr.es %*% yus[,ii]
        proposed.log.likelihood = proposed.log.likelihood - 0.5 * t(yus[,ii]) %*% proposed.inv.Corr.es %*% yus[,ii]
      }
      log.mh.ratio = proposed.log.likelihood - current.log.likelihood	
      if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
      if(log(runif(1))<log.mh.ratio)
      {
        be[tt] = current.be[tt] = proposed.be[tt]
        Corr.es = current.Corr.es = proposed.Corr.es
        inv.Corr.es = current.inv.Corr.es = proposed.inv.Corr.es
      }	
    }
    
    for(tt in 1:((d-2)*(d-1)/2))
    {
      proposed.thetae = current.thetae
      proposed.thetae[tt] = sample(c(current.thetae[tt],current.thetae[tt]+c(-1,1)*2*3.14/40),size=1)
      while(abs(proposed.thetae[tt])>3.14)
        proposed.thetae[tt] = sample(c(current.thetae[tt],current.thetae[tt]+c(-1,1)*2*0.99/40),size=1)		
      proposed.Corr.es = formCorr.es(be,proposed.thetae)
      proposed.inv.Corr.es = solve(proposed.Corr.es)
      
      proposed.log.likelihood = 0		
      current.log.likelihood = 0
      for(ii in 1:N)
      {
        current.log.likelihood = current.log.likelihood - 0.5 * t(yus[,ii]) %*% current.inv.Corr.es %*% yus[,ii]
        proposed.log.likelihood = proposed.log.likelihood - 0.5 * t(yus[,ii]) %*% proposed.inv.Corr.es %*% yus[,ii]
      }
      log.mh.ratio = proposed.log.likelihood - current.log.likelihood	
      if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
      if(log(runif(1))<log.mh.ratio)
      {
        thetae[tt] = current.thetae[tt] = proposed.thetae[tt]
        Corr.es = current.Corr.es = proposed.Corr.es
        inv.Corr.es = current.inv.Corr.es = proposed.inv.Corr.es
      }	
    }
    
    
    
    
    if((iii>burnin/2)&&(T))
    {
      ### Updating xs (and us)
      
      if(q>0)
      {
        for(jj in 1:q)
        {
          proposed.xs[jj,] = rtnorm(n,mean=current.xs[jj,],sd=0.1,lower=x.lwr[jj],upper=x.upr[jj])	
          TempMat = abs(matrix(rep(proposed.xs[jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[jj,],n),n,x.grid.length,byrow=T))
          close.ind = apply(TempMat,1,which.min) 
          
          proposed.prior[jj,] = B.basis.store[jj,close.ind,]%*%B.basis.density.coeffs(thetas.x.density[jj,],delta.t[jj])
          current.prior[jj,] = d.ordinates.x[jj,]
          
          proposed.xs.eps[jj,] = B.basis.store[jj,close.ind,]%*%thetas.eps[jj,]
          proposed.xs.consumption.days[jj,] = proposed.xs[jj,]/pnorm(proposed.xs.eps[jj,])
          proposed.xs.consumption.days[jj,proposed.xs.consumption.days[jj,]>max(knots.t[jj,])] = max(knots.t[jj,])-range.start.xs[jj]/1000  # For fixing numerical issues
          
          TempMat = abs(matrix(rep(proposed.xs.consumption.days[jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[jj,],n),n,x.grid.length,byrow=T))
          close.ind = apply(TempMat,1,which.min) 
          proposed.vars[jj,] = B.basis.store[jj,close.ind,]%*%exp(thetas[jj,])
          proposed.vars[jj,proposed.vars[jj,]==0] = max(proposed.vars[jj,])
          
          proposal.p.to.c[jj,] = dtnorm(current.xs[jj,],mean=proposed.xs[jj,],sd=diff(range(start.xs[jj,]))/6,lower=x.lwr[jj],upper=x.upr[jj])
          proposal.c.to.p[jj,] = dtnorm(proposed.xs[jj,],mean=current.xs[jj,],sd=diff(range(start.xs[jj,]))/6,lower=x.lwr[jj],upper=x.upr[jj])
          
          temp.xsrs = Fx.norm.B.basis(xs[jj,],knots.t[jj,],thetas.x.density[jj,])
          temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
          temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
          yxs[jj,] = qnorm(temp.xsrs)
          
          proposed.yxs.temp = Fx.norm.B.basis(proposed.xs[jj,],knots.t[jj,],thetas.x.density[jj,])
          proposed.yxs.temp[proposed.yxs.temp<0.001] = 0.001
          proposed.yxs.temp[proposed.yxs.temp>0.999] = 0.999
          proposed.yxs[jj,] = qnorm(proposed.yxs.temp)
          
          proposed.us[jj,] = ws[jj,]-proposed.xs.consumption.days[jj,inds]
          
          proposed.us.eps[jj,] = ws.eps[jj,]-proposed.xs.eps[jj,inds]		
          temp.current.marginal.likelihood.eps = dnorm(current.us.eps[jj,],mean=0,sd=1)
          temp.proposed.marginal.likelihood.eps = dnorm(proposed.us.eps[jj,],mean=0,sd=1)
          current.marginal.likelihood.eps[jj,] = tapply(temp.current.marginal.likelihood.eps,inds,"prod")
          proposed.marginal.likelihood.eps[jj,] = tapply(temp.proposed.marginal.likelihood.eps,inds,"prod")
        }
      }
      
      if(p>0)
      {
        for(jj in 1:p)
        {
          proposed.xs[q+jj,] = rtnorm(n,mean=current.xs[q+jj,],sd=diff(range(start.xs[q+jj,]))/6,lower=x.lwr[q+jj],upper=x.upr[q+jj])
          TempMat = abs(matrix(rep(proposed.xs[q+jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[q+jj,],n),n,x.grid.length,byrow=T))
          close.ind = apply(TempMat,1,which.min) 
          
          k.x = max(z.x[jj,])		
          proposed.prior[q+jj,] = fx_mixtnorm(proposed.xs[q+jj,],pi.x[jj,1:k.x],mu.x[1:k.x],sqrt(sigmasq.x[1:k.x]),x.lwr[q+jj],x.upr[q+jj]) 
          current.prior[q+jj,] = fx_mixtnorm(current.xs[q+jj,],pi.x[jj,1:k.x],mu.x[1:k.x],sqrt(sigmasq.x[1:k.x]),x.lwr[q+jj],x.upr[q+jj])
          
          proposed.xs.consumption.days[q+jj,] = proposed.xs[q+jj,]
          
          TempMat = abs(matrix(rep(proposed.xs[q+jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[q+jj,],n),n,x.grid.length,byrow=T))
          close.ind = apply(TempMat,1,which.min) 
          proposed.vars[q+jj,] = B.basis.store[q+jj,close.ind,]%*%exp(thetas[q+jj,])
          
          proposal.p.to.c[q+jj,] = dtnorm(current.xs[q+jj,],mean=proposed.xs[q+jj,],sd=diff(range(start.xs[q+jj,]))/6,lower=x.lwr[q+jj],upper=x.upr[q+jj])
          proposal.c.to.p[q+jj,] = dtnorm(proposed.xs[q+jj,],mean=current.xs[q+jj,],sd=diff(range(start.xs[q+jj,]))/6,lower=x.lwr[q+jj],upper=x.upr[q+jj])
          
          proposed.yxs.temp = Fx_mixtnorm(proposed.xs[q+jj,],pi.x[jj,1:k.x],mu.x[1:k.x],sqrt(sigmasq.x[1:k.x]),x.lwr[q+jj],x.upr[q+jj])
          proposed.yxs.temp[proposed.yxs.temp<0.001] = 0.001
          proposed.yxs.temp[proposed.yxs.temp>0.999] = 0.999
          proposed.yxs[q+jj,] = qnorm(proposed.yxs.temp)
          
          proposed.us[q+jj,] = ws[q+jj,]-proposed.xs.consumption.days[q+jj,inds]
        }
      }
      
      current.copula.prior = as.vector(dmvnorm(t(current.yxs),rep(0,d),Corr.xs) * exp(0.5*colSums(current.yxs^2)))
      proposed.copula.prior = as.vector(dmvnorm(t(proposed.yxs),rep(0,d),Corr.xs)  * exp(0.5*colSums(proposed.yxs^2)))
      
      for(jj in 1:d)
      {
        k.u = max(z.u[jj,])	
        temp.current.marginal.likelihood = fu_mixnorm(current.us[jj,],mean=0,sd=sqrt(current.vars[jj,inds]),pi.u[jj,1:k.u],params.u[1:k.u,])
        temp.proposed.marginal.likelihood = fu_mixnorm(proposed.us[jj,],mean=0,sd=sqrt(proposed.vars[jj,inds]),pi.u[jj,1:k.u],params.u[1:k.u,])
        current.marginal.likelihood[jj,] = tapply(temp.current.marginal.likelihood,inds,"prod")
        proposed.marginal.likelihood[jj,] = tapply(temp.proposed.marginal.likelihood,inds,"prod")
        
        proposed.yus.temp = Fu_mixnorm(proposed.us[jj,],0,sqrt(proposed.vars[jj,inds]),pi.u[jj,1:k.u],params.u[1:k.u,])
        proposed.yus.temp[proposed.yus.temp<0.001] = 0.001
        proposed.yus.temp[proposed.yus.temp>0.999] = 0.999
        proposed.yus[jj,] = qnorm(proposed.yus.temp)	
      }
      
      if(q>0)
      {
        for(jj in 1:q)
        {
          current.marginal.likelihood[jj,] = current.marginal.likelihood.eps[jj,] * current.marginal.likelihood[jj,] 		
          proposed.marginal.likelihood[jj,] = proposed.marginal.likelihood.eps[jj,] * proposed.marginal.likelihood[jj,]
        }
      }
      
      temp.current.copula.likelihood = dmvnorm(t(current.yus),rep(0,d),Corr.es) * exp(0.5*colSums(current.yus^2))
      temp.proposed.copula.likelihood = dmvnorm(t(proposed.yus),rep(0,d),Corr.es) * exp(0.5*colSums(proposed.yus^2))
      current.copula.likelihood = as.vector(tapply(temp.current.copula.likelihood,inds,"prod"))
      proposed.copula.likelihood = as.vector(tapply(temp.proposed.copula.likelihood,inds,"prod"))
      
      mh.ratio = (apply(proposal.p.to.c,2,'prod') * apply(proposed.prior,2,'prod') * proposed.copula.prior * apply(proposed.marginal.likelihood,2,'prod') * proposed.copula.likelihood)/(apply(proposal.c.to.p,2,'prod') * apply(current.prior,2,'prod') * current.copula.prior * apply(current.marginal.likelihood,2,'prod') * current.copula.likelihood)
      
      mh.ratio[is.nan(mh.ratio)] = 0
      
      u = runif(n)
      inds.to.replace = (1:n)[u<mh.ratio]
      xs[,inds.to.replace] = current.xs[,inds.to.replace] = proposed.xs[,inds.to.replace]
      vars[,inds.to.replace] = current.vars[,inds.to.replace] = proposed.vars[,inds.to.replace]
      us = current.us = ws - xs.consumption.days[,inds]
      
      if(q>0){xs.eps[,inds.to.replace] = current.xs.eps[,inds.to.replace] = proposed.xs.eps[,inds.to.replace]
      xs.consumption.days[,inds.to.replace] = current.xs.consumption.days[,inds.to.replace] = proposed.xs.consumption.days[,inds.to.replace]
      us.eps = current.us.eps = ws.eps-xs.eps[,inds]}
      
      
      
      for(jj in 1:d)
      {
        ### Updating theta
        
        proposed.thetas[jj,] = t(rmvnorm(1,current.thetas[jj,],diag(rep(sig.tune.thetas.1,(K.t+1)))+sig.tune.thetas.2*prop.sig.thetas[jj,,]))
        proposed.thetas[jj,1] = 2*proposed.thetas[jj,2]-proposed.thetas[jj,3]
        TempMat = abs(matrix(rep(xs.consumption.days[jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[jj,],n),n,x.grid.length,byrow=T))
        close.ind = apply(TempMat,1,which.min) 
        proposed.vars[jj,] = B.basis.store[jj,close.ind,]%*%exp(proposed.thetas[jj,])
        
        current.log.prior = - t(current.thetas[jj,])%*%P.t%*%current.thetas[jj,]/(2*s2t[jj])
        proposed.log.prior = - t(proposed.thetas[jj,])%*%P.t%*%proposed.thetas[jj,]/(2*s2t[jj])
        
        temp.current.likelihood = d.restricted.mix.norm(us[jj,],mean=rep(0,times=N),sd=sqrt(current.vars[jj,inds]),params.u[z.u[jj,],])
        temp.proposed.likelihood = d.restricted.mix.norm(us[jj,],mean=rep(0,times=N),sd=sqrt(proposed.vars[jj,inds]),params.u[z.u[jj,],])
        
        temp.current.likelihood[temp.current.likelihood==0] = 10^(-100)  # Numerical Stability
        temp.proposed.likelihood[temp.proposed.likelihood==0] = 10^(-100)  # Numerical Stability
        current.log.likelihood = sum(log(temp.current.likelihood))
        proposed.log.likelihood = sum(log(temp.proposed.likelihood))
        
        log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
        if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
        if(log(runif(1))<log.mh.ratio)
        {
          thetas[jj,] = current.thetas[jj,] = proposed.thetas[jj,]
          vars[jj,] = current.vars[jj,] = proposed.vars[jj,]
        }
        
        ### Updating s2t
        
        s2t[jj] = 1/rgamma(1,shape=alpha.t+(K.t+1)/2,rate=beta.t+t(thetas[jj,])%*%P.t%*%thetas[jj,])
      }
    }
    
    
    
    if(p>0)	
    {       
      ### Updating xs.trans (to be used in next iteration)
      
      for(jj in 1:p)
      {
        xs.trans.temp.num = pnorm((xs[q+jj,]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))-pnorm((x.lwr[q+jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))
        xs.trans.temp.den = pnorm((x.upr[q+jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))-pnorm((x.lwr[q+jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))
        xs.trans.temp.den[xs.trans.temp.den<=0] = 0.001  # Numerical Stability
        xs.trans.temp = xs.trans.temp.num/xs.trans.temp.den
        xs.trans.temp[xs.trans.temp>=1] = 0.999  # Numerical Stability
        xs.trans.temp[xs.trans.temp<=0] = 0.001  # Numerical Stability
        xs.trans.temp = mu.x[z.x[jj,]]+sqrt(sigmasq.x[z.x[jj,]])*qnorm(xs.trans.temp)
        xs.trans.temp[xs.trans.temp<(x.lwr[jj]-10)] = x.lwr[jj] - 10  # Numerical Stability
        xs.trans.temp[xs.trans.temp>(x.upr[jj]+10)] = x.upr[jj] + 10  # Numerical Stability
        xs.trans[jj,] = xs.trans.temp
      }
    }
    
    
    
    if(q>0)	
    {
      for(jj in 1:q)
      {
        ### Updating ws.eps
        
        xs.star = xs.eps[jj,inds]
        indices_not_consumed_temp = which(indices_not_consumed[jj,]==1)	
        ws.eps[jj,indices_not_consumed_temp] = rtnorm(num_proxies_not_consumed[jj],mean=xs.star[indices_not_consumed_temp],sd=1,lower=-Inf,upper=0)
        indices_consumed_temp = which(indices_consumed[jj,]==1)	
        ws.eps[jj,indices_consumed_temp] = rtnorm(num_proxies_consumed[jj],mean=xs.star[indices_consumed_temp],sd=1,lower=0,upper=Inf)   
        
        ### Updating thetas.eps
        
        TempMat = abs(matrix(rep(xs[jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[jj,],n),n,x.grid.length,byrow=T))
        close.ind = apply(TempMat,1,which.min) 
        Mu.thetas.eps = numeric(K.t+1)
        Sigma.thetas.eps = matrix(0,K.t+1,K.t+1)
        for(ii in 1:n)
        {
          for(hh in inds1[ii]:inds2[ii])
            Mu.thetas.eps = Mu.thetas.eps + ws.eps[jj,hh] * B.basis.store[jj,close.ind[ii],]
          Sigma.thetas.eps = Sigma.thetas.eps + mis[ii] * B.basis.store[jj,close.ind[ii],] %*% t(B.basis.store[jj,close.ind[ii],])
        }
        Sigma.thetas.eps = solve(P.t/s2t.eps[jj]+solve(Sigma0.thetas.eps)+Sigma.thetas.eps)
        Mu.thetas.eps = Sigma.thetas.eps %*% (solve(Sigma0.thetas.eps)%*%Mu0.thetas.eps[jj,] + Mu.thetas.eps)
        thetas.eps[jj,] = t(rmvnorm(1,Mu.thetas.eps,Sigma.thetas.eps))
        
        xs.eps[jj,] = current.xs.eps[jj,] = B.basis.store[jj,close.ind,]%*%thetas.eps[jj,]
        xs.consumption.days[jj,] = current.xs.consumption.days[jj,] = xs[jj,]/pnorm(xs.eps[jj,])
        xs.consumption.days[jj,xs.consumption.days[jj,]>max(knots.t[jj,])] = max(knots.t[jj,])-range.start.xs[jj]/1000  # For fixing numerical issues
        
        TempMat = abs(matrix(rep(xs.consumption.days[jj,],x.grid.length),n,x.grid.length)-matrix(rep(x.grid[jj,],n),n,x.grid.length,byrow=T))
        close.ind = apply(TempMat,1,which.min) 
        vars[jj,] = current.vars[jj,] = B.basis.store[jj,close.ind,]%*%exp(thetas[jj,])
        
        ### Updating s2t.eps
        
        s2t.eps[jj] = 1/rgamma(1,shape=alpha.t.eps[jj]+(K.t+1)/2,rate=beta.t.eps[jj]+t(thetas.eps[jj,])%*%P.t%*%thetas.eps[jj,])
        
        ### Updating ws when not observed
        
        C2.u = rbinom(N,1,prob=c(params.u[z.u[jj,],1],1-params.u[z.u[jj,],1]))+1
        mu.u = cbind(params.u[,2]*(1-params.u[,1])/sqrt(params.u[,1]^2+(1-params.u[,1])^2),-params.u[,2]*params.u[,1]/sqrt(params.u[,1]^2+(1-params.u[,1])^2))	
        mu.w = xs.consumption.days[jj,inds] + sqrt(vars[jj,inds])*mu.u[cbind(z.u[jj,],C2.u)]
        sigmasq.w = vars[jj,inds]*params.u[cbind(z.u[jj,],C2.u+2)]
        indices_not_consumed_temp = which(indices_not_consumed[jj,]==1)	
        ws[jj,indices_not_consumed_temp] = rnorm(num_proxies_not_consumed[jj],mean=mu.w[indices_not_consumed_temp],sd=sqrt(sigmasq.w[indices_not_consumed_temp]))
      }
      
      
      ### Re-updating us.eps and us
      
      us.eps = current.us.eps = ws.eps-xs.eps[,inds]
      us = current.us = (ws-xs.consumption.days[,inds])
    }
    
    
    ### Save results after burn-in
    
    if(iii>burnin)
    {
      if(q>0)
      {
        for(jj in 1:q)
          density.x.est[jj,] = density.x.est[jj,] + B.basis.store[jj,,]%*%B.basis.density.coeffs(thetas.x.density[jj,],delta.t[jj])
      }
      if(p>0)
      {
        for(jj in 1:p)
        {
          k.x = max(z.x[jj,])
          for(kk in 1:k.x)
            density.x.est[q+jj,] = density.x.est[q+jj,] + pi.x[jj,kk]*dtnorm(x.grid[q+jj,],mu.x[kk],sqrt(sigmasq.x[kk]),x.lwr[q+jj],x.upr[q+jj])
          if(sum(is.nan(density.x.est))>0)
            stop("NAN occurred!")
        }
      }
      
      for(jj in 1:d)
      {
        k.u = max(z.u[jj,])
        var.e[jj] = var.e.fn(pi.u[jj,1:k.u],params.u[1:k.u,])
        density.e.est[jj,] = density.e.est[jj,] + d.scaled.restricted.mix.norm(e.grid,0,1,pi.u[jj,1:k.u],params.u[1:k.u,])
        var.est[jj,] = var.est[jj,] + B.basis.var.grid.knots.t[jj,,] %*% exp(thetas[jj,]) * var.e[jj]
      }
      
      if(q>0)
      {
        for(jj in 1:q)
        {
          prob.consumption.est[jj,] = prob.consumption.est[jj,] + pnorm(B.basis.store[jj,,]%*%thetas.eps[jj,])
        }
      }
      Corr.xs.est = Corr.xs.est + Corr.xs	
      Corr.es.est = Corr.es.est + Corr.es			
    }
  }
  
  
  
  #################################  
  ### Compute Density Estimates ###
  #################################
  
  Corr.xs.est = Corr.xs.est/(simsize-burnin)
  Corr.xs.est
  Corr.es.est = Corr.es.est/(simsize-burnin)
  Corr.es.est
  if(q>0)
    prob.consumption.est = prob.consumption.est/(simsize-burnin)
  for(jj in 1:d)
  {
    density.x.est[jj,] = density.x.est[jj,]/(simsize-burnin)
    var.true[jj,] = abs(var.grid[jj,]/3)^2
    var.est[jj,] = var.est[jj,]/(simsize-burnin)
  }
  for(jj in 1:d)
    density.e.est[jj,] = density.e.est[jj,]/(simsize-burnin)
  
  s2is.new = wbars.new = matrix(0,d,n)
  for(jj in 1:d)
  {
    wbars.new[jj,] = tapply(ws[jj,],inds,"mean")  
    s2is.new[jj,] = as.vector(tapply(ws[jj,],inds,var))
  }
  
  
  
  ####################################  
  ### Compute 2D Density Estimates ###
  ####################################
  
  CDF.x.est <- y.est.grid <- matrix(0,nrow=d,ncol=x.grid.length)
  for(jj in 1:d)
  {
    CDF.x.est[jj,] <- cumsum(density.x.est[jj,]) * (x.grid[jj,2]-x.grid[jj,1])
    CDF.x.est[jj,CDF.x.est[jj,]>1] <- 1
    y.est.grid[jj,] <- qnorm(CDF.x.est[jj,])
  }
  y.est.grid[y.est.grid==Inf] <- 3.65
  y.est.grid[y.est.grid==-Inf] <- -3.65
  
  indices <- combn(d,2)
  y.est.grid.2D <- array(0,dim=c(2,x.grid.length*x.grid.length))
  density.x.est.2D <- array(0,dim=c(2,x.grid.length*x.grid.length))
  density.x.est.2D.array <- array(0,dim=c(d,d,x.grid.length,x.grid.length))
  for(kk in 1:ncol(indices))
  {
    dim.indices <- indices[,kk]
    d1 <- dim.indices[1]
    d2 <- dim.indices[2]
    
    density.x.est.2D[1,] <- rep(density.x.est[d1,],times=x.grid.length)
    density.x.est.2D[2,] <- rep(density.x.est[d2,],each=x.grid.length)
    
    y.est.grid.2D <- rbind(rep(y.est.grid[d1,],times=x.grid.length),rep(y.est.grid[d2,],each=x.grid.length))
    
    density.x.est.2D.array.vect <- (dmvnorm(t(y.est.grid.2D),rep(0,2),Corr.xs.est[c(d1,d2),c(d1,d2)]) / dmvnorm(t(y.est.grid.2D),rep(0,2),diag(1,2))) * apply(density.x.est.2D,2,prod)
    density.x.est.2D.array.vect[is.na(density.x.est.2D.array.vect)] <- 0
    density.x.est.2D.array[d1,d2,,] <- matrix(density.x.est.2D.array.vect,nrow=x.grid.length,byrow=T)
  }
  
  
  
  
  
  ###################################################
  ### Compute the Densities of Normalized Intakes ###
  ###################################################
  
  CDF.x.est <- y.est.grid <- matrix(0,nrow=d,ncol=x.grid.length)
  for(jj in 1:d)
  {
    CDF.x.est[jj,] <- cumsum(density.x.est[jj,]) * (x.grid[jj,2]-x.grid[jj,1])
    CDF.x.est[jj,CDF.x.est[jj,]>1] <- 1
    y.est.grid[jj,] <- qnorm(CDF.x.est[jj,])
  }
  xs.nrmlzd.est <- matrix(1,d,n)
  x.nrmlzd.grid <- matrix(1,d,x.grid.length)
  for(jj in 1:(d-1))
  {
    xs.nrmlzd.est[jj,] <- xs[jj,]/xs[d,]
    x.nrmlzd.grid[jj,] <- seq(0,quantile(xs.nrmlzd.est[jj,],0.9999),len=x.grid.length)
  }
  density.x.nrmlzd.est <- matrix(0,d-1,x.grid.length)
  y.nrmlzd.est <- matrix(0,d,x.grid.length)
  y.nrmlzd.est[d,] <- y.est.grid[d,]
  density.x.nrmlzd.est <- density.x.nrmlzd.est.tmp1 <- density.x.est
  for(jj in 1:(d-1))
  {
    for(ll in 1:x.grid.length)
    {
      x.jj.grid.ll <- x.nrmlzd.grid[jj,ll]*x.grid[d,]
      TempMat <- abs(matrix(rep(x.jj.grid.ll,x.grid.length),x.grid.length,x.grid.length)-matrix(rep(x.grid[jj,],n),x.grid.length,x.grid.length,byrow=T))
      close.ind <- apply(TempMat,1,which.min)
      
      y.nrmlzd.est[jj,] <- y.est.grid[jj,close.ind]
      density.x.nrmlzd.est.tmp1[jj,] <- density.x.est[jj,close.ind] 
      density.x2.nrmlzd.est.tmp2 <- (dmvnorm(t(y.nrmlzd.est[c(jj,d),]),rep(0,2),Corr.xs.est[c(jj,d),c(jj,d)]) / dmvnorm(t(y.nrmlzd.est[c(jj,d),]),rep(0,2),diag(1,2))) * apply(density.x.nrmlzd.est.tmp1[c(jj,d),],2,prod)
      density.x2.nrmlzd.est.tmp2[which(is.nan(density.x2.nrmlzd.est.tmp2))] <- 0     
      density.x.nrmlzd.est[jj,ll] <- sum(x.grid[d,]*density.x2.nrmlzd.est.tmp2) * (x.grid[jj,2]-x.grid[jj,1])
    }
  }
  
  
  
  
  ####################  
  ### Plot Results ###
  ####################
  
  if(Plot_Results==TRUE)
  {
    filename = paste("HTR","d",d,sep="_")
    for(jj in 1:d)
      filename = paste(filename,var.names[jj],sep="_")
    filename = paste(filename,"simsize",simsize,"burnin",burnin,"seedno",seedno,sep="_")
    filename = paste(filename,".pdf",sep="")
    print(filename)
    pdf(filename)
    
    
    
    ### Plots Similar to Figures 5 & 8 in the Main Paper ###
    
    par(mar=rep(2,4),mfrow=c(d,3))
    layout(matrix(c(1,1,2,1,1,3,4,4,5,4,4,6,7,7,8,7,7,9,10,10,11,10,10,12),nrow=d*2,ncol=3,byrow=T))
    for(jj in 1:d) 
    {
      # Plot estimated marginal densities of variable of interest x
      plot(x.grid[jj,],density.x.est[jj,],xlab="x",ylab=expression(f[x]),xlim=c(0,x.upr[jj]),ylim=c(0,max(density.x.est[jj,])),type="l",lwd=1)
      
      # Plot estimated marginal densities of the scaled errors e
      plot(e.grid,density.e.est[jj,],xlab="e",ylab=expression(f[e]),ylim=c(0,max(density.e.est[jj,])),type="l",lwd=1)
      
      # Plot estimated variance functions superimposed over subject-specific sample means-vs-variances
      plot(wbars.new[jj,],s2is.new[jj,],xlab="x",ylab="v(x)",xlim=c(0,quantile(wbars.new[jj,],0.98)),ylim=c(0,quantile(s2is.new[jj,],0.98)),pch="*",col="gray50")
      points(var.grid[jj,1:70],var.est[jj,1:70],type="l",lwd=1,lty=1,col="black")
    }
    par(mfrow=c(1,1))
    
    
    
    ### Plots Similar to Figures 6 & 9 in the Main Paper ###
    
    if(q>0)
    {
      par(mar=rep(2,4),mfrow=c(q,1))
      for(jj in 1:q) 
      {
        # Plot estimated probabilties of reporting a positive recall as a function of x
        max.index = x.grid.length-length(which(x.grid[jj,]-min(max(x.grid[jj,]),max(knots.t[jj,]))>0))
        plot(x.grid[jj,1:max.index],prob.consumption.est[jj,1:max.index],ylim=c(0,1),type="l",lwd=1,xlab="x",ylab="prob",lty=1)
        title(main="Probability of consumption")
      }
    }
    
    
    
    ### Plots Similar to Figures S.4 in the Supplementary Materials ###
    
    laymat <- diag(1:d)
    laymat[lower.tri(laymat)] <-  (d+1):(d*(d+1)/2)
    laymat <- t(laymat) 										# fixes the upper part
    laymat_temp <- diag(rep(0,d))
    laymat_temp[upper.tri(laymat_temp)] <- ((d*(d+1)/2)+1):d^2
    laymat <- laymat+t(laymat_temp)								# fixes the lower part
    laymat
    layout(laymat) 
    
    par(mar=rep(1.85,4))
    my.cols <- rev(brewer.pal(11, "RdYlBu"))
    x.grid.lims <- c(1:100*x.grid.length/100)
    for(jj in 1:d) 
    {
      # Plot estimated marginal densities of variable of interest x in the diagonal panels
      plot(x.grid[jj,x.grid.lims],density.x.est[jj,x.grid.lims],ylim=c(0,max(density.x.est[jj,x.grid.lims])),col="blue",type="l",lwd=2,lty=1)
    }
    for(i in 1:(d-1))
      for(j in (i+1):d)
      {
        # Plot estimated bivariate joint densities of variable of interest x in the off-diagonal panels
        image(x.grid[i,x.grid.lims],x.grid[j,x.grid.lims],density.x.est.2D.array[i,j,x.grid.lims,x.grid.lims],xlim=c(0,x.upr[i]),ylim=c(0,x.upr[j]),col=terrain.colors(100))
        contour(x.grid[i,x.grid.lims],x.grid[j,x.grid.lims],density.x.est.2D.array[i,j,x.grid.lims,x.grid.lims],nlevels=45,drawlabels=FALSE,col=my.cols[6:11],add=T)
        legend("bottomright",legend=paste(i,j),cex=1,bty="n")
      }
    
    
    
    ### Plots Similar to Figures 7 & 10 in the Main Paper ###

    par(mar=rep(2,4),mfrow=c(d-1,1))
    for(jj in 1:(d-1)) 
    {
      xlwrjj = quantile(xs.nrmlzd.est[jj,],0.001)
      xuprjj = quantile(xs.nrmlzd.est[jj,],0.999)
      ylimjj = max(density.x.nrmlzd.est[jj,])
      # Plot estimated marginal densities of normalized intakes
      plot(x.nrmlzd.grid[jj,],density.x.nrmlzd.est[jj,],xlim=c(xlwrjj,xuprjj),ylim=c(0,ylimjj),type="l",lwd=1,main="")
      title(main="Density of normalized intakes")
    } 
    par(mfrow=c(1,1))
    
    dev.off()
  }
  
  
  
  
  
  #############################  
  ### Save Plotting Results ###
  #############################
  
  filename = paste("HTR","d",d,sep="_")
  for(jj in 1:d)
    filename = paste(filename,var.names[jj],sep="_")
  filename = paste(filename,"simsize",simsize,"burnin",burnin,"seedno",seedno,sep="_")
  filename = paste(filename,".RData",sep="")
  print(filename)
  if(Save_Results==TRUE)
    save(file=filename,xs,knots.t,x.grid,density.x.est,density.x.est.2D.array,e.grid,density.e.est,var.grid,var.est,prob.consumption.est,xs.nrmlzd.est,x.nrmlzd.grid,density.x.nrmlzd.est)
  }
  
  
  
  
  