require(FME);
################################## SubROUTINE used in this script#######################################################
# The "cost" function 
FRAUZcost <- function (pars,y0,SS,times,obs) {
	pars <- pars*SS; 
	out <- ode(y=y0,parms=pars,times=times,method='ode45',func='derivs',initfunc='initmod',dllname='flauz',nout=4);
	cost_1 <- modCost(model=out, obs=obs, weight='std')
	return(cost_1)
}
#model 
FRAUZ<-function(y0,pars,SS,times){
	pars=pars*SS;
	out <- ode(y=y0,parms=pars,times=times,method='ode45',func='derivs',initfunc='initmod',dllname='flauz',nout=4);
	return(out);
}
mdMCMC<-function(ff=FRAUZcost,pp,obs,niter=35000,y=y0,SS,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
	pp <- pp*SS;
	Fit <- modFit(f=ff,p=pp,y=y0,SS=SS,times=times,obs=obs,lower=lower,upper=upper);
	pp <- Fit$par;
	
	MCMC1 <- modMCMC(f=ff, p=pp, niter=10000, y=y0,SS,times=times,obs=obs,
	 wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
	
	Cov0<-cov(MCMC1$pars)*2.4^2/(length(pp));
	lower=(summary(MCMC1)[5,])*0.75;
	upper=summary(MCMC1)[7,]*1.5;
	pp=unlist(summary(MCMC1)[6,]);
	MCMC <- modMCMC(f=ff, p=pp, niter=niter, y=y0,SS,times=times,obs=obs,jump=Cov0,
	 wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper, burninlength=burninlength);
	 
	 return(MCMC);
}

FRAUZcost1 <- function (pars,y0,times,obs) {
	out <- ode(y=y0,parms=pars,times=times,method='ode45',func='derivs',initfunc='initmod',dllname='flauz',nout=4);
	cost_1 <- modCost(model=out, obs=obs, weight='std')
	return(cost_1)
}
#model 
FRAUZ1<-function(y0,pars,times){
	out <- ode(y=y0,parms=pars,times=times,method='ode45',func='derivs',initfunc='initmod',dllname='flauz',nout=4);
	return(out);
}
#plot

mdMCMC1<-function(ff=FRAUZcost1,pp,obs,niter=35000,y=y0,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
	Fit <- modFit(f=ff,p=pp,y=y0,times=times,obs=obs,lower=lower,upper=upper);
	pp <- Fit$par;
	MCMC1 <- modMCMC(f=ff, p=pp, niter=10000, y=y0,times=times,obs=obs,wvar0=wvar0, updatecov=updatecov,
		lower=lower, upper=upper,burninlength=1000);
	
	Cov0<-cov(MCMC1$pars)*2.4^2/(length(pp));
	lower=(summary(MCMC1)[5,])*0.75;
	upper=summary(MCMC1)[7,]*1.5;
	pp=unlist(summary(MCMC1)[6,]);
	MCMC <- modMCMC(f=ff, p=pp, niter=niter, y=y0,times=times,obs=obs,jump=Cov0,wvar0=wvar0, updatecov=updatecov,
		lower=lower, upper=upper, burninlength=burninlength);
	 
	 return(MCMC);
}

valid_plotMCMCpdf<-function(obs,mod_a,fn='valid_plot.pdf'){
	tn<-names(mod_a);
	nob<-ncol(obs);
	yl<-apply(mod_a,FUN=range,2)
	nm<-ncol(mod_a)-2;
	n<-floor(sqrt(nm-1));
	m<-floor((nm-1)/(n-1));
	times<-obs[,'time'];
	par0<-par();
	pdf(fn);
		par(mfcol=c(4,5), mar=c(2,3,2,0.5))
		for (i in 2:(length(tn)-2)){
			plot(x=times,y=mod_a[,i],type='l',lwd=2,xlab='',ylab='',ylim=c(yl[1,i]*0.9,yl[2,i]*1.1));title(tn[i]);
			if(i<=(nob)) points(times,obs[,i],type='p',cex=2,pch=19);
		};	
	par(par0);
	dev.off();
	return(1);
}
valid_plotMCMC<-function(obs,mod_a){
	tn<-names(mod_a);
	nob<-ncol(obs);
	yl<-apply(mod_a,FUN=range,2)
	nm<-ncol(mod_a)-2;
	n<-floor(sqrt(nm-1));
	m<-floor((nm-1)/(n-1));
	times<-obs[,'time'];
	par0<-par();
	win.graph();
		par(mfcol=c(4,5), mar=c(2,3,2,0.5))
		for (i in 2:(length(tn)-2)){
			plot(x=times,y=mod_a[,i],type='l',lwd=2,xlab='',ylab='',ylim=c(yl[1,i]*0.9,yl[2,i]*1.1));title(tn[i]);
			if(i<=(nob)) points(times,obs[,i],type='p',cex=2,pch=19);
		};	
	par(par0);
	return(1);
}
