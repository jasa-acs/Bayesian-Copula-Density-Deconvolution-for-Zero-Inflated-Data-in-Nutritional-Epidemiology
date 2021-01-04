
#####################################
### Bayesian Copula Deconvolution ###
#####################################

# Codes accompanying "Bayesian Copula Density Deconvolution for Zero-Inflated Data with Applications in Nutritional Epidemiology" by Sarkar, Pati, Mallick and Carroll.
# Codes written by Abhra Sarkar (abhra.sarkar@utexas.edu), last modified on Dec 15, 2019, in Austin, TX

# The current file is for univariate density deconvolution for variables with strictly continuously measured surrogates. 
# The method uses mixtures of truncated-normals with shared atoms to model the density of interest, 
# mixtures of moment-restricted normals to model the density of the measurement errors,
# and mixtures of B-spines to model the conditional variability of the measurement errors.
# See paper for additional details.



#############
### Input ###
#############

# While running from within the file 'Bayes_Copula_Decon_MVT.R' that implements the multivariate method, these arguments are read from the original file.
# The univariate method can also be independently implemented using the current file.  

# ws <- strictly continuously measured surrogate values
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

# knots <- knot-points for constructing the B-splines bases that model the conditional variability of the measurement errors 
# thetas <- estimated coefficients of B-splines bases that model the conditional variability of the measurement errors 
# xs <- estimated subject-specific values of the variable of interest
# us <-  estimated subject and replicate-specific values of the measurement errors
# z.xs <- mixture component labels for the mixture model for the density of interest
# pi.xs <- mixture component probabilities for the mixture model for the density of interest
# params.xs <- mixture component-specific parameters for the mixture model for the density of interest
# z.us <- mixture component labels for the mixture model for the density of the measurement errors
# pi.us <- mixture component probabilities for the mixture model for the density of the measurement errors
# params.us <- mixture component-specific parameters for the mixture model for the density of the measurement errors



UNIV_DECON_REGULAR = function(ws, xs.lwr, xs.upr, mis, z.xs.max, z.us.max, K.t, simsize, burnin, show_progress=TRUE, plot_results=TRUE)
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
  current.xs = start.xs = xs
  us = ws - rep(xs,times=mis)
  range.start.xs = diff(range(xs))
  s2is = as.vector(tapply(ws,inds,var))
  xs.grid = seq(xs.lwr,xs.upr,length=500)
  xs.grid.length = length(xs.grid)
  
  alpha.xs = 1
  
  # Normal 
  mu0.xs = mean(xs)
  sigmasq0.xs = var(xs)
  
  
  # Inverse-Gamma (Independnt from mu - independence is important)
  a.sigmasq.xs = 1
  b.sigmasq.xs = 1
  
  pi.xs = rep(1/z.xs.max,z.xs.max)
  clusters.xs = kmeans(xs,z.xs.max) 
  mu.xs = clusters.xs$center
  z.xs = clusters.xs$cluster
  sigmasq.xs = rep(var(xs)/5,z.xs.max)
  
  d.ordinates.xs = matrix(0,nrow=n,ncol=z.xs.max)
  
  ### Prior and initialization of s2t and thetas
  alpha.t = 100
  beta.t = 1
  s2t = 0.01
  P.t = P.mat(K.t+1) 	# penalty matrix
  knots.t = seq(xs.lwr,xs.upr,length=K.t)
  
  optim_results = optim(rep(1,K.t+1), fr, NULL, xs, mis, knots.t, P.t, s2t, us, method = "BFGS")
  thetas = current.thetas = start.thetas = optim_results$par
  prop.sig.thetas = make.positive.definite(prop.sig.thetas.fn(thetas,xs,mis,us,s2t,K.t,P.t,knots.t,n))
  
  var.grid = seq(xs.lwr,xs.upr,length=100)
  vars = current.vars = t(B.basis(xs,knots.t)%*%exp(current.thetas))
  B.basis.var.grid.knots.t = B.basis(var.grid,knots.t)
  B.basis.store = B.basis(xs.grid,knots.t)
  close.ind = rep(0,n)
  
  ### Prior and initialization for mixture
  simsize.mh.us = 10
  z.us = rep(1,N)
  alpha.us = 0.1
  params.us = matrix(c(0.5,0,1,1),nrow=z.us.max,ncol=4,byrow=T)		# unique values
  pi.us = rep(1/z.us.max,z.us.max)
  d.ordinates.us = matrix(0,nrow=N,ncol=z.us.max)
  
  
  
  	#########################
  	### Tuning Parameters ###
  	#########################
  
  
  sig.tune.thetas.1 = 0
  sig.tune.thetas.2 = 0.1
  
  
  
  
  
  	###############################
  	### Storage for MCMC Output ###
  	###############################
  
  
  es.grid = seq(-3,3,length=500)
  density.xs.est = numeric(xs.grid.length)
  var.es = numeric(1)
  var.est = numeric(length(var.grid))
  density.es.est = numeric(length(es.grid))
  prob.consumption.est = numeric(xs.grid.length)
  
  proposed.xs = current.xs = xs 
  proposed.us = current.us = us
  proposed.vars = current.vars = vars 
  
  current.likelihood = proposed.likelihood = matrix(1,2,n)
  temp.proposed.us.likelihood = temp.current.us.likelihood = matrix(1,2,N)
  
  thetas.est  = numeric(length(thetas))
  thetas.MCMC = matrix(0,nrow=simsize,ncol=length(thetas))
  
  
  
  	##################
  	### Start MCMC ###
  	##################
  
  
  for (iii in 1:simsize)
  	{
  	if((show_progress==TRUE)&&(iii%%10==0))
  		print(iii)
  
  
  	### Updating z.xs
  	for(kk in 1:z.xs.max)
  		d.ordinates.xs[,kk] = dtnorm(xs,mu.xs[kk],sqrt(sigmasq.xs[kk]),lower=xs.lwr,upper=xs.upr)
  	d.ordinates.xs[is.nan(d.ordinates.xs)] = 0
  	d.ordinates.xs[is.infinite(d.ordinates.xs)] = 0
  	for(ii in 1:n)
  		z.xs[ii] = sample(z.xs.max,1,prob=pi.xs*d.ordinates.xs[ii,])
  
  	
  	### Updating cluster probabilities
  	n.kk.xs = tabulate(z.xs,nbins=z.xs.max)
  	pi.xs = rdirichlet(1,alpha.xs/z.xs.max+n.kk.xs)	
  
  
  	### Updating mu.xs, sigmasq.xs
  	
  	xs.trans = mu.xs[z.xs]+sqrt(sigmasq.xs[z.xs])*qnorm((pnorm((xs-mu.xs[z.xs])/sqrt(sigmasq.xs[z.xs]))-pnorm((xs.lwr-mu.xs[z.xs])/sqrt(sigmasq.xs[z.xs])))/(pnorm((xs.upr-mu.xs[z.xs])/sqrt(sigmasq.xs[z.xs]))-pnorm((xs.lwr-mu.xs[z.xs])/sqrt(sigmasq.xs[z.xs]))))
  	xs.trans[xs.trans < xs.lwr - 10] = xs.lwr - 10
  	xs.trans[xs.trans > xs.upr + 10] = xs.upr + 10
  
  	for(kk in 1:z.xs.max)
  		{
  		temp = which(z.xs==kk)
  		xspool = xs.trans[temp]
  
  		sigmasq.temp = 1/(n.kk.xs[kk]/sigmasq.xs[kk] + 1/sigmasq0.xs)
  		mu.temp = (sum(xspool)/sigmasq.xs[kk] + mu0.xs/sigmasq0.xs) * sigmasq.temp
  		mu.xs[kk] = rnorm(1,mu.temp,sqrt(sigmasq.temp))	
  						
  		post.a.sigmasq.xs = a.sigmasq.xs + length(xspool)/2
  		post.b.sigmasq.xs = b.sigmasq.xs + sum((xspool-mu.xs[kk])^2)/2
  		sigmasq.xs[kk] = 1/rgamma(1,shape=post.a.sigmasq.xs,rate=post.b.sigmasq.xs)
  		}
  
  
  
  	### Updating xs (and us)
  
  	proposed.xs = rtnorm(n,mean=current.xs,sd=0.1,lower=xs.lwr,upper=xs.upr)	
    TempMat = abs(matrix(rep(proposed.xs,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
  	close.ind = apply(TempMat,1,which.min) 
    proposed.vars = B.basis.store[close.ind,]%*%exp(thetas)
  		
  	proposed.prior = dtnorm(proposed.xs,mu.xs[z.xs],sqrt(sigmasq.xs[z.xs]),lower=xs.lwr,upper=xs.upr)
  	current.prior = dtnorm(current.xs,mu.xs[z.xs],sqrt(sigmasq.xs[z.xs]),lower=xs.lwr,upper=xs.upr)
  
  	proposed.us = ws-rep(proposed.xs,times=mis)
  		
  	k.us = max(z.us)	
  	temp.current.us.likelihood = fu_mixnorm(current.us,mean=0,sd=rep(sqrt(current.vars),times=mis),pi.us[1:k.us],params.us[1:k.us,])
  	temp.proposed.us.likelihood = fu_mixnorm(proposed.us,mean=0,sd=rep(sqrt(proposed.vars),times=mis),pi.us[1:k.us],params.us[1:k.us,])
  	current.likelihood = tapply(temp.current.us.likelihood,inds,"prod")
  	proposed.likelihood = tapply(temp.proposed.us.likelihood,inds,"prod")
  
  	mh.ratio = (proposed.prior * proposed.likelihood * dtnorm(current.xs,mean=proposed.xs,sd=0.1,lower=xs.lwr,upper=xs.upr)) / (current.prior * current.likelihood * dtnorm(proposed.xs,mean=current.xs,sd=0.1,lower=xs.lwr,upper=xs.upr))
  
  	mh.ratio[is.nan(mh.ratio)] = 0
  
  	u = runif(n)
  	inds.to.replace = (1:n)[u<mh.ratio]
  	xs[inds.to.replace] = current.xs[inds.to.replace] = proposed.xs[inds.to.replace]
  	vars[inds.to.replace] = current.vars[inds.to.replace] = proposed.vars[inds.to.replace]
  
  	us = current.us = ws - rep(xs,times=mis)
  
  
  
  	### Updating thetas
  
  	proposed.thetas = t(rmvnorm(1,current.thetas,(diag(rep(sig.tune.thetas.1,(K.t+1)))+sig.tune.thetas.2*prop.sig.thetas)))
    TempMat = abs(matrix(rep(xs,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
  	close.ind = apply(TempMat,1,which.min) 
    proposed.vars = B.basis.store[close.ind,]%*%exp(proposed.thetas)
  		
  	current.log.prior = - t(current.thetas)%*%P.t%*%current.thetas/(2*s2t)
  	proposed.log.prior = - t(proposed.thetas)%*%P.t%*%proposed.thetas/(2*s2t)
  
  	temp.current.likelihood = d.restricted.mix.norm(us,mean=rep(0,times=N),sd=rep(sqrt(current.vars),times=mis),params.us[z.us,])
  	temp.proposed.likelihood = d.restricted.mix.norm(us,mean=rep(0,times=N),sd=rep(sqrt(proposed.vars),times=mis),params.us[z.us,])
  
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
  	
  	
  	### Updating z.us
  		
  	for(ii in 1:N)
  		{			
  		prob.us = pi.us * d.restricted.mix.norm(us[ii],mean=0,sd=sqrt(vars[inds[ii]]), params.us)
  		if(sum(prob.us)==0) {prob.us=rep(1/z.us.max,z.us.max)}
  		z.us[ii] = sample(1:z.us.max,1,TRUE,prob.us)   # New z.us[ii] drawn
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
  			temp.proposed.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),proposed.params.us))
  			temp.current.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),params.us[kk,]))
  			temp.proposed.log.likelihood[is.infinite(temp.proposed.log.likelihood)] = 0
  			temp.current.log.likelihood[is.infinite(temp.current.log.likelihood)] = 0
  			proposed.log.likelihood = sum(temp.proposed.log.likelihood)
  			current.log.likelihood = sum(temp.current.log.likelihood)
  
  			log.acc.prob = proposed.log.likelihood-current.log.likelihood
  			if(log(runif(1))<log.acc.prob)
  				params.us[kk,] = proposed.params.us
  			}
  		}
  	if(k.us<z.us.max)
  	for(kk in (k.us+1):z.us.max)
  		params.us[kk,] = r.proposal.params.restricted.mix.norm(1,1,3,3,3,3,3)	
  
  	var.es = var.e.fn(pi.us[1:k.us],params.us[1:k.us,])
  	params.us[1:k.us,2] = params.us[1:k.us,2]/sqrt(var.es)
  	params.us[1:k.us,3:4] = params.us[1:k.us,3:4]/var.es
  	
  
  
  	thetas.MCMC[iii,] = thetas
  
  	if(iii>burnin)
  		{		
  		for(kk in 1:z.xs.max)
  			density.xs.est = density.xs.est + pi.xs[kk]*dtnorm(xs.grid,mu.xs[kk],sqrt(sigmasq.xs[kk]),lower=xs.lwr,upper=xs.upr)
  		
  		k.us = max(z.us)
  		var.es = var.e.fn(pi.us[1:k.us],params.us[1:k.us,])
  		density.es.est = density.es.est + d.scaled.restricted.mix.norm(es.grid,0,1,pi.us[1:k.us],params.us[1:k.us,])
  		var.est = var.est + B.basis.var.grid.knots.t %*% exp(thetas) * var.es
  		thetas.est = thetas.est + log(var.es) + thetas
  	  }
  	
  	}
  
  
  
  density.xs.est = density.xs.est/(simsize-burnin)
  density.es.est = density.es.est/(simsize-burnin)
  var.est = var.est/(simsize-burnin)
  thetas.est = thetas.est/(simsize-burnin)
  
  thetas.final = thetas.est
  xs.final = xs
  var.final = sqrt(B.basis(xs.final,knots.t)%*%exp(thetas.final))
  us.final = (ws-rep(xs.final,times=mis))
  
  
  if(plot_results==TRUE)
  	{
  	dev.new()
  	par(mfrow=c(2,2))
  	plot(xs.grid,density.xs.est,xlab="x",ylab="f(x)",type="l",lty=1,col="green3",lwd=3)
  
  	plot(es.grid,density.es.est,xlab="e",ylab="f(e)",type="l",lty=1,col="green3",lwd=3)
  	points(es.grid,dnorm(es.grid),type="l",lty=1)
  
  	plot(xs,s2is,pch="*",xlab="x",ylab="v(x)")
  	points(var.grid,var.est,type="l",lty=1,col="blue",lwd=2)
  	points(var.grid,B.basis.var.grid.knots.t%*%exp(thetas.final),type="l",lty=1,col="green3",lwd=2)
  	par(mfrow=c(1,1))
    }
  
  
  params.xs = cbind(mu.xs,sigmasq.xs)
  return(list(knots=knots.t, thetas=thetas.final, 
  	xs=xs.final, us=us.final, 
  	z.xs=z.xs, pi.xs=pi.xs, params.xs=params.xs, 
  	z.us=z.us, pi.us=pi.us, params.us=params.us))
}

