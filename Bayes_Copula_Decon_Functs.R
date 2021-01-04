#############################################################


newlog <- function(x) { (log(x) > 0)*log(x) + (log(x) < 0)*((x -1) - (x-1)^2/2 + (x-1)^3/3  - (x-1)^4/4) }


#############################################################


P.mat <- function(K)
	{
	# penalty matrix for density
	D <- diag(rep(1,K))
	D <- diff(diff(D))
	P <- t(D)%*%D 
	return(P)
	}


#############################################################


B.basis <- function(x,knots)
	{
	delta <- knots[2]-knots[1]
	n <- length(x)
	K <- length(knots)
	B <- matrix(0,n,K+1)
	for (jj in 1:(K-1))
      		{
		act.inds <- (1:n)[(x>=knots[jj])&(x<=knots[jj+1])]
		act.x <- x[act.inds]
		resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
        
		B[act.inds,jj] <- (1/2)*(1-resc.x)^2
		B[act.inds,jj+1] <- -(resc.x^2)+resc.x+1/2
		B[act.inds,jj+2] <- (resc.x^2)/2
		}
	return(B)
	}


#############################################################


B.basis.normalized <- function(x,knots)		# equidistant knots
	{
	B <- B.basis(x,knots)
	K <- length(knots)
	delta <- knots[2]-knots[1]
	
	B[,1] <- B[,1] /(delta/6)
	B[,2] <- B[,2] /(5*delta/6)
	B[,K] <- B[,K] /(5*delta/6)
	B[,K+1] <- B[,K+1] /(delta/6)
	B[,3:(K-1)] <- B[,3:(K-1)]/delta
	
	return(B)
	}


#############################################################


B.basis.local.integral <- function(x,knots)
  {
	delta <- knots[2]-knots[1]
	n <- length(x)
	K <- length(knots)
	B <- matrix(0,n,K+1)
	for (jj in 1:(K-1))
    {
		act.inds <- (1:n)[(x>=knots[jj])&(x<=knots[jj+1])]
		act.x <- x[act.inds]
		resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
        
		B[act.inds,jj] <- 5*delta/6 + (delta/2)*(resc.x-resc.x^2+resc.x^3/3)
		B[act.inds,jj+1] <- delta/6 + delta*(-resc.x^3/3+resc.x^2/2+resc.x/2)
		B[act.inds,jj+2] <- delta*resc.x^3/6
		}
	act.inds <- (1:n)[(x>=knots[1])&(x<=knots[2])]
	B[act.inds,1] <- B[act.inds,1] - 5*delta/6 
	B[act.inds,2] <- B[act.inds,2] - delta/6 
	act.inds <- (1:n)[(x>=knots[2])&(x<=knots[3])]
	B[act.inds,2] <- B[act.inds,2] - delta/6 
	return(B)
  }


B.basis.cum.integral <- function(x,knots)
  {
	delta <- knots[2]-knots[1]
	n <- length(x)
	K <- length(knots)
	B <- matrix(0,n,K+1)
	for (jj in 1:(K-1))
    {
		act.inds <- (1:n)[(x>=knots[jj])]
		act.x <- x[act.inds]
		resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
    resc.x[resc.x>1] <- 1
    
		B[act.inds,jj] <- 5*delta/6 + (delta/2)*(resc.x-resc.x^2+resc.x^3/3)
		B[act.inds,jj+1] <- delta/6 + delta*(-resc.x^3/3+resc.x^2/2+resc.x/2)
		B[act.inds,jj+2] <- delta*resc.x^3/6
	  }
	act.inds <- (1:n)[(x>=knots[1])]
	B[act.inds,1] <- B[act.inds,1] - 5*delta/6 
	B[act.inds,2] <- B[act.inds,2] - delta/6 
	return(B)
  }


#############################################################


B.basis.density.coeffs <- function(thetas.density,delta)
  {
  coeffs <- exp(thetas.density)/((delta/6) * (exp(thetas.density[1])+5*exp(thetas.density[2])+6*sum(exp(thetas.density[3:(K.t-1)]))+5*exp(thetas.density[K.t])+exp(thetas.density[K.t+1])))
  return(coeffs)
  }


#############################################################


Fx.norm.B.basis <- function(x,knots,thetas.density)
  {
	delta <- knots[2]-knots[1]
	Fx <- B.basis.cum.integral(x,knots)%*%B.basis.density.coeffs(thetas.density,delta)
	return(Fx)
  }


#############################################################


rnormbspline <- function(n,coeffs,minx,maxx)
	{
	grid <- seq(minx,maxx,length=500)
	delta <- (grid[2]-grid[1])
	knots <- seq(minx,maxx,len=(length(coeffs)-1))
	y <- B.basis.normalized(grid,knots)%*%(coeffs*length(knots))
	y <- y/(sum(y)*delta)
	Fy <- cumsum(y)*delta
	u <- runif(n)	
	x <- numeric(n)
	for(ii in 1:n)
		x[ii] <- grid[min(which(u[ii]<Fy))+1]
	x <- x + runif(n,0,delta)
	return(x)
	}
	

#############################################################


dnormbspline <- function(x,coeffs,minx,maxx)	# x must be equidistant grid
	{
	knots <- seq(minx,maxx,len=(length(coeffs)-1))
	y <- B.basis.normalized(x,knots)%*%(coeffs*length(knots))
	delta <- (x[2]-x[1])
	y <- y/(sum(y)*delta)
	return(y)
	}
	

#############################################################


rskewnorm <- function(n,mean,sd,skewness)
	{
	c = 2/3.1415926
	delta = skewness/sqrt(1+skewness^2)
	xi = - (sqrt(c) * skewness)/sqrt(1 + skewness^2 * (1-c))
	omega = sqrt(1+xi^2)

	y = delta * abs(rnorm(n,0,1)) + (1-delta^2)^(1/2) * rnorm(n,0,1)
	y = xi + omega * y
	return(mean+sd*y)
	}


dskewnorm <- function(x,mean,sd,skewness)
	{
	c = 2/3.1415926
	delta = skewness/sqrt(1+skewness^2)
	zeta1 = delta * sqrt(c)
	zeta2 = sqrt(1 - c*delta^2)
	y = numeric(length(x))
	xmod = zeta1 + zeta2*(x-mean)/sd
	for(i in 1:length(x))
      	y[i] = (2*zeta2/sd[i]) * dnorm(xmod[i]) * pnorm(skewness*xmod[i])
	return(y)
	}


#############################################################


d.scaled.restricted.mix.norm <- function(x,mean,sd,pi,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)

  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	sd.e = sqrt(var.e.fn(pi,params))
	y = matrix(0,nrow=length(p),ncol=length(x))
	for(kk in 1:length(p))
	    	y[kk,] = pi[kk]*(p[kk]*dnorm(x,(mean+sd*mu1[kk])/sd.e,sd*sqrt(sigmasq1[kk])/sd.e) + (1-p[kk])*dnorm(x,(mean+sd*mu2[kk])/sd.e,sd*sqrt(sigmasq2[kk])/sd.e))
	y = colSums(y)
	return(y)
	}



d.restricted.mix.norm <- function(x,mean,sd,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

  y = p*dnorm(x,mean+sd*mu1,sd*sqrt(sigmasq1)) + (1-p)*dnorm(x,mean+sd*mu2,sd*sqrt(sigmasq2))
	return(y)
	}



mixnorm <- function(n,pi,p,mu_curl,sigmasq1,sigmasq2,e.grid,plot=TRUE)
	{
	m = length(pi)
	y = numeric(n)
	density <- numeric(length(e.grid))
	inds = sample(1:m,n,TRUE,prob=pi)
	for(ii in 1:m)
		{
		temp = which(inds==ii)
		y[temp] = r.restricted.mix.norm(length(temp),c(p[ii],mu_curl[ii],sigmasq1[ii],sigmasq2[ii]))
		}
	params = cbind(p,mu_curl,sigmasq1,sigmasq2)
	y = y/sqrt(var.e.fn(pi,params))
	density <- d.scaled.restricted.mix.norm(e.grid,mean=0,sd=1,pi,params)
	if(plot)
		{
		hist(y,xlim=c(min(e.grid),max(e.grid)),breaks=30,freq=FALSE)
		points(e.grid,density,type="l")
		}

	return(list(es=y,density=density))
	}


r.restricted.mix.norm <- function(n,params)
	{
	p = params[1]
	mu_curl = params[2]
	sigmasq1 = params[3]
	sigmasq2 = params[4]

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)

  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

  inds = sample(0:1,n,TRUE,prob=c(p,(1-p)))
	y = numeric(n)
	temp = which(inds==0)
	y[temp] = rnorm(length(temp),mu1,sqrt(sigmasq1)) 
	temp = which(inds==1)
	y[temp] = rnorm(length(temp),mu2,sqrt(sigmasq2))
	return(y)
	}



var.e.fn <- function(pi,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)


  mu1 = c1*mu_curl
	mu2 = c2*mu_curl


  y = p*(mu1^2+sigmasq1) + (1-p)*(mu2^2+sigmasq2)
	y = sum(pi*y)
	return(y)
	}



gamma.e.fn <- function(pi,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)


  mu1 = c1*mu_curl
	mu2 = c2*mu_curl


  m2 = p*(mu1^{2}+sigmasq1) + (1-p)*(mu2^{2}+sigmasq2)
	m2 = sum(pi*m2)

  m3 = p*(mu1^{3}+3*mu1*sigmasq1) + (1-p)*(mu2^{3}+3*mu2*sigmasq2)
	m3 = sum(pi*m3)

  m4 = p*(mu1^{4}+6*mu1^{2}*sigmasq1+3*sigmasq1^{2}) + (1-p)*(mu2^{3}+6*mu2^{2}*sigmasq2+3*sigmasq2^{2})
	m4 = sum(pi*m4)

	gamma1 = m3/m2^{3/2}
	gamma2 = m4/m2^{2}-3

	return(c(gamma1,gamma2))
	}


r.proposal.params.restricted.mix.norm <- function(p.a,p.b,sigmasq.mu_curl,s11,s12,s21,s22)
	{
	p = rbeta(1,p.a,p.b)

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
	mu_curl = rnorm(1,0,sqrt(sigmasq.mu_curl))
	sigmasq1 = 1/rgamma(1,s11,s12)
	sigmasq2 = 1/rgamma(1,s21,s22)

	y = c(p,mu_curl,sigmasq1,sigmasq2)
	return(y)
	}    


r.tnorm.proposal.params.restricted.mix.norm <- function(params)
	{
	current.p <- params[1]
	current.mu_curl <- params[2]
	current.sigmasq1 <- params[3]
	current.sigmasq2 <- params[4]

	p = rtnorm(1,current.p,0.01,lower=0,upper=1)
	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
    
	sigmasq1 = rtnorm(1,current.sigmasq1,0.1,lower=0,upper=Inf)
	sigmasq2 = rtnorm(1,current.sigmasq2,0.1,lower=0,upper=Inf)

  mu_curl = rnorm(1,current.mu_curl,0.1)

	y = c(p,mu_curl,sigmasq1,sigmasq2)
	return(y)
	}    




#############################################################

dlaplace <- function(x,mu,b)
	{
	y = exp(-abs(x-mu)/b)/(2*b)
	return(y)
	}

plaplace <- function(x,mu,b)
	{
	y = numeric(length(x))
	for(ii in 1:length(y))
		{
		if(x[ii]<mu) 
			y[ii] = (1/2)*exp(-abs(x[ii]-mu)/b)
		else
			y[ii] = 1-(1/2)*exp(-abs(x[ii]-mu)/b)
		}		
	return(y)
	}


rlaplace <- function(n,mu,b)
	{
	u = runif(n)
	y = mu - b*sign(u-0.5)*log(1-2*abs(u-0.5))
	return(y)
	}



mixlaplace <- function(n,pi,mu,b,e.grid,plot=TRUE)		# Produces scaled mixtures of Laplace
	{
	m = length(pi)
	y = numeric(n)
	density = numeric(length(e.grid))
	inds = sample(1:m,n,TRUE,pi)
    mu.e = sum(pi*mu)
	sd.e = sqrt(sum(pi*(2*b^{2}+mu^{2}))-mu.e^{2})
	for(ii in 1:m)
		{
		temp = which(inds==ii)
		y[temp] = rlaplace(length(temp),mu[ii],b[ii])
		density = density + pi[ii]*dlaplace(e.grid,(mu[ii]-mu.e)/sd.e,b[ii]/sd.e)
		}
	y = (y-mu.e)/sd.e
	if(plot)
		{
		hist(y,xlim=c(min(e.grid),max(e.grid)),freq=FALSE,breaks=30)
		points(e.grid,density,type="l")
		}
	return(list(es=y,density=density))
	}


d.scaled.restricted.mixlaplace <- function(x,pi,mu,b)		# Produces scaled mixtures of Laplace
	{
	y = matrix(0,nrow=length(pi),ncol=length(x))
    mu.e = sum(pi*mu)
	sd.e = sqrt(sum(pi*(2*b^{2}+mu^{2}))-mu.e^{2})
	for(kk in 1:length(pi))
	    y[kk,] = pi[kk]*dlaplace(x,(mu[kk]-mu.e)/sd.e,b[kk]/sd.e)
	y = colSums(y)
	return(y)
	}



gamma.e.Laplace.fn <- function(pi,mu,b)
	{
	m1.dash = sum(pi*mu) 
	m2.dash = sum(pi*(2*b^{2}+mu^{2}))
	m3.dash = sum(pi*(mu^{3}+6*b^{2}*mu))
	m4.dash = sum(pi*(mu^{4}+12*b^{2}*mu^{2}+24*b^{4}))

	m2 = m2.dash-m1.dash^{2}
  m3 = m3.dash - 3*m2.dash*m1.dash + 2*m1.dash^{3}
	m4 = m4.dash - 4*m3.dash*m1.dash + 6*m2.dash*m1.dash^{2} - 3*m1.dash^{4}

  gamma1 = m3/m2^{3/2}
	gamma2 = m4/m2^{2}-3

	return(c(gamma1,gamma2))
	}



#############################################################


fr <- function(thetas,xs,mis,knots.t,P.t,s2t,us)			# function to be maximized
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)
	y = - t(thetas) %*% P.t %*% thetas / (2*s2t) - sum(rep(log(vars),times=mis))/2 - sum(us^2/rep(vars,times=mis))/2
	return(-y)
  }

fr.density <- function(thetas,xs,knots.t,P.t,s2t,density.est)			# function to be maximized
	{
  delta = knots.t[2]-knots.t[1]
	density.est2 = B.basis(xs,knots.t)%*%B.basis.density.coeffs(thetas,delta)
	y = - t(thetas) %*% P.t %*% thetas / (2*s2t) - sum((density.est-density.est2)^2)
	return(-y)
	}

fr.prob.consumption <- function(thetas,xs,knots.t,P.t,s2t,probs.emp)			# function to be maximized
{
  probs.est = as.vector(pnorm(B.basis(xs,knots.t)%*%thetas))
  y = - t(thetas) %*% P.t %*% thetas / (2*s2t) - sum((probs.emp-probs.est)^2)
  return(-y)
}


############################################################


gr <- function(thetas,xs,mis,knots.t,P.t,s2t,us)			# gradient function of fr
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)
	y = - P.t %*% thetas / s2t
	B.basis.components = matrix(0,nrow=K.t+1,ncol=n)
	for(kk in 1:(K.t+1))
		{
   		thetas.new = rep(0,K.t+1)
		thetas.new[kk] = 1
		B.basis.components[kk,] = B.basis(xs,knots.t)%*%thetas.new
		}
	for(kk in 1:(K.t+1))
		for(ii in 1:n)
			for(jj in (sum(mis[1:ii-1])+1):sum(mis[1:ii]))
				y[kk] = y[kk] - (1 - (us[jj]^2)/vars[ii]) * B.basis.components[kk,ii]*exp(thetas[kk])/(2*vars[ii])
	return(-y)
	}


############################################################


prop.sig.thetas.fn <- function(thetas,xs,mis,us,s2t,K.t,P.t,knots.t,n)
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)

	B.basis.components = matrix(0,nrow=K.t+1,ncol=n);
	for(kk in 1:(K.t+1))
		{
   	thetas.new = rep(0,K.t+1)
		thetas.new[kk] = 1
		B.basis.components[kk,] = B.basis(xs,knots.t)%*%thetas.new
	  }

	prop.sig.thetas = matrix(0,nrow=K.t+1,ncol=K.t+1)
	for(kk in 1:(K.t+1))
		{
		for(ll in kk:(K.t+1))
			{
			if(kk==ll)
			  {
	  			for(ii in 1:n)
	     				for(jj in (sum(mis[1:(ii-1)])+1):sum(mis[1:ii]))
	        				prop.sig.thetas[kk,ll] = prop.sig.thetas[kk,ll] + ((us[jj]^2)/vars[ii] - 1/2)*(B.basis.components[kk,ii]*B.basis.components[ll,ii])*exp(thetas[kk]+thetas[ll])/(vars[ii]^2) + (1-(us[jj]^2)/vars[ii])*B.basis.components[kk,ii]*exp(thetas[kk])/(2*vars[ii])
			  }
			else
			  {
	  			for(ii in 1:n)
	     				for(jj in (sum(mis[1:(ii-1)])+1):sum(mis[1:ii]))
	       					prop.sig.thetas[kk,ll] = prop.sig.thetas[kk,ll] + ((us[jj]^2)/vars[ii] - 1/2)*(B.basis.components[kk,ii]*B.basis.components[ll,ii])*exp(thetas[kk]+thetas[ll])/(vars[ii]^2)
			  }
		  } 
		}
	for(kk in 2:(K.t+1))
   		for(ll in 1:(kk-1))
			    prop.sig.thetas[kk,ll] = prop.sig.thetas[ll,kk]
	prop.sig.thetas = prop.sig.thetas + P.t/s2t
	prop.sig.thetas = round(solve(prop.sig.thetas),4)
	
	return(prop.sig.thetas)
	}


#############################################################


pi.fn <- function(x,alpha,beta,xstar,K)
	{
	w <- pnorm(alpha-beta*abs(x-xstar)^2)
	pi <- numeric(K)
	pi[1] <- w[1]
	for(k in 2:(K-1))
		pi[k] <- w[k]*prod(1-w[1:(k-1)])
	pi[K] <- prod(w[1:(K-1)])
	return(pi)
	}


#############################################################
# Folded Normal Distribution 

dfldnorm <- function(x,mu,sd) 
	{
	dens = ifelse(x<0,0,dnorm(-x,mu,sd) + dnorm(x,mu,sd))
	return(dens)
	}



#############################################################


fx_mixnorm <- function(x,w,m,s) 
	{
	y = matrix(0,nrow=length(w),ncol=length(x))
	for(kk in 1:length(w))
	    y[kk,] = w[kk]*dnorm(x,m[kk],s[kk])
	y = colSums(y)
	return(y)
	}
	
Fx_mixnorm <- function(x,w,m,s) 
	{
	y = matrix(0,nrow=length(w),ncol=length(x))
	for(kk in 1:length(w))
	    y[kk,] = w[kk]*pnorm(x,m[kk],s[kk])
	y = colSums(y)
	return(y)
	}
	
Fx_mixnorm_inv <- function(p,w,m,s,br=c(-1000,1000))	# Currently accepts only scalar p
	{
	G <- function(x) {Fx_mixnorm(x,w,m,s) - p}
	return(uniroot(G,br)$root) 
	}
	
	
	
fx_mixtnorm <- function(x,w,m,s,lwr,upr) 
	{
	y = matrix(0,nrow=length(w),ncol=length(x))
	for(kk in 1:length(w))
	    y[kk,] = w[kk]*dtnorm(x,m[kk],s[kk],lwr,upr)
	y = colSums(y)
	return(y)
	}
	
Fx_mixtnorm <- function(x,w,m,s,lwr,upr) 
	{
	y = matrix(0,nrow=length(w),ncol=length(x))
	for(kk in 1:length(w))
	    y[kk,] = w[kk]*ptnorm(x,m[kk],s[kk],lwr,upr)
	y = colSums(y)
	return(y)
	}
	
Fx_mixtnorm_inv <- function(p,w,m,s,lwr,upr,br=c(-1000,1000))	# Currently accepts only scalar p
	{
	G <- function(x) {Fx_mixtnorm(x,w,m,s,lwr,upr) - p}
	return(uniroot(G,br)$root) 
	}
	
	

Fe_scaled_mixnorm <- function(x,pi,params) 
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	sd.e = sqrt(var.e.fn(pi,params))
	y = matrix(0,nrow=length(p),ncol=length(x))
	for(kk in 1:length(p))
	    y[kk,] = pi[kk]*(p[kk]*pnorm(x,mu1[kk]/sd.e,sqrt(sigmasq1[kk])/sd.e) + (1-p[kk])*pnorm(x,mu2[kk]/sd.e,sqrt(sigmasq2[kk])/sd.e))
	y = colSums(y)
	return(y)
	}
		
Fe_scaled_mixnorm_inv <- function(p,pi,params,br=c(-100,100))	# Currently accepts only scalar p
	{
	G <- function(x) {Fe_scaled_mixnorm(x,pi,params) - p}
	return(uniroot(G,br)$root) 
	}


Fe_scaled_laplace <- function(x,pi,mu,b) 
	{
	y = matrix(0,nrow=length(pi),ncol=length(x))
  mu.e = sum(pi*mu)
	sd.e = sqrt(sum(pi*(2*b^{2}+mu^{2}))-mu.e^{2})
	for(kk in 1:length(pi))
	    y[kk,] = pi[kk]*plaplace(x,(mu[kk]-mu.e)/sd.e,b[kk]/sd.e)
	y = colSums(y)
	return(y)
	}
		
Fe_scaled_laplace_inv <- function(p,pi,mu,b,br=c(-100,100))
	{
	G <- function(x) {Fe_scaled_laplace(x,pi,mu,b) - p}
	return(uniroot(G,br)$root) 
	}


fu_mixnorm <- function(x,mean,sd,pi,params) 
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	y = matrix(0,nrow=length(p),ncol=length(x))
	for(kk in 1:length(p))
	  y[kk,] = pi[kk]*(p[kk]*dnorm(x,mean+sd*mu1[kk],sd*sqrt(sigmasq1[kk])) + (1-p[kk])*dnorm(x,mean+sd*mu2[kk],sd*sqrt(sigmasq2[kk])))
	y = colSums(y)
	return(y)
	}


Fu_mixnorm <- function(x,mean,sd,pi,params) 
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	y = matrix(0,nrow=length(p),ncol=length(x))
	for(kk in 1:length(p))
	    y[kk,] = pi[kk]*(p[kk]*pnorm(x,mean+sd*mu1[kk],sd*sqrt(sigmasq1[kk])) + (1-p[kk])*pnorm(x,mean+sd*mu2[kk],sd*sqrt(sigmasq2[kk])))
	y = colSums(y)
	return(y)
	}


fu_eps_mixnorm <- function(x,mean,sd,pi,params) 
	{
	if(is.matrix(params))
		{
		mu = params[,1]
		sigmasq = params[,2]
		}
	else
		{
		mu = params[1]
		sigmasq = params[2]
		}
    
	y = matrix(0,nrow=length(pi),ncol=length(x))
	for(kk in 1:length(pi))
	    y[kk,] = pi[kk]*dnorm(x,mean+sd*mu[kk],sd*sqrt(sigmasq[kk])) 
	y = colSums(y)
	return(y)
	}


Fu_eps_mixnorm <- function(x,mean,sd,pi,params) 
	{
	if(is.matrix(params))
		{
		mu = params[,1]
		sigmasq = params[,2]
		}
	else
		{
		mu = params[1]
		sigmasq = params[2]
		}
    
	y = matrix(0,nrow=length(pi),ncol=length(x))
	for(kk in 1:length(pi))
	    y[kk,] = pi[kk]*pnorm(x,mean+sd*mu[kk],sd*sqrt(sigmasq[kk]))
	y = colSums(y)
	return(y)
	}



#############################################################


formCorr.xs <- function(b,theta)
	{
	p = length(b)+1	
	V = matrix(0,p,p)
	V[1,1] = 1
	V[2,1] = b[1]
	V[2,2] = sqrt(1-b[1]*b[1])
	for(jj in 3:p)
		{
		q1jj = 1+(jj-3)*(jj-2)/2
		q2jj = (jj-2)*(jj-1)/2
		V[jj,1] = b[jj-1]*sin(theta[q1jj])
		if(jj>3)
			for(kk in 2:(jj-2))
				V[jj,kk] = b[jj-1]*prod(cos(theta[q1jj:(q1jj+kk-2)]))*sin(theta[q1jj+kk-1])
		V[jj,jj-1] = b[jj-1]*prod(cos(theta[q1jj:q2jj]))	
		V[jj,jj] = sqrt(1-b[jj-1]*b[jj-1])	
		}
	R = V%*%t(V)
	return(R)	
	}
	

formCorr.es <- function(b,theta)
	{
	p = length(b)+1	
	V = matrix(0,p,p)
	V[1,1] = 1
	V[2,1] = b[1]
	V[2,2] = sqrt(1-b[1]*b[1])
	for(jj in 3:p)
		{
		q1jj = 1+(jj-3)*(jj-2)/2
		q2jj = (jj-2)*(jj-1)/2
		V[jj,1] = b[jj-1]*sin(theta[q1jj])
		if(jj>3)
			for(kk in 2:(jj-2))
				V[jj,kk] = b[jj-1]*prod(cos(theta[q1jj:(q1jj+kk-2)]))*sin(theta[q1jj+kk-1])
		V[jj,jj-1] = b[jj-1]*prod(cos(theta[q1jj:q2jj]))	
		V[jj,jj] = sqrt(1-b[jj-1]*b[jj-1])	
		}
	R = V%*%t(V)
	return(R)	
	}


