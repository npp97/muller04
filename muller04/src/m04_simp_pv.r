library("FME")
library('compiler')

dyn.load('flauz.dll')
loadcmp('flauz_subs.Rc')

#######################################################################################################################
#------------------Inti STATE VARIABLES
y0 = c(NH4 = 44.45,
	   NO3 = 25.43,
	   NH4_N15 = 4.48,
	   NO3_N15 = 0.033,
	   TON_N15=0.049,
	   TON = 37.0,	
	   Mib = 0,
	   PR = 17,
	   Hum = 20,
	   Mib_N15 = 0.049,
	   Hum_N15 = 0.049,
	   NH3 = 0.000,
	   N2O = 0.000,
	   NH3_N15 = 0.000,
	   N2O_N15 = 0.000,
	   PR_N15 = 0
		   );
#--------------Parameters
pars0 = c(C_hu_nh4 = 0.00001, 
		C_micb_hum= 0., 
		C_pr_nh4= 0. , 
		C_pr_micb= 0.0001,
		C_micb_nh4=0.0001,
		K_nh4mic = 0.02,
		K_no3mic= 0.0001,
		K_nh4no3= 0.04,
		K_no3ng= 0.0001,
		K_nh4nh3= 0.0001,
		K_no3nh4=0
		);
#----------Switch to on/off some processes,the values must be 0 or 1.
SS =  c(S_hu_nh4 = 1, 
		S_micb_hum= 1, 
		S_pr_nh4= 1, 
		S_pr_micb= 1,
		S_micb_nh4=1,
		S_nh4mic = 1,
		S_no3mic= 1,
		S_nh4no3= 1,
		S_no3ng= 1,
		S_nh4nh3= 1,
		S_no3nh4= 1
		);


#Read data into memory, variables must be some of the following variables:
#-----------------------------------------------------------------------------------------------------------------------------------------------
#Time, NH4, Err_NH4, NO3, Err_NO3, NH4_N15, Err_NH4_N15, NO3_N15, Err_NO3_N15, ON, Err_ON, ON_N15, Err_ON_N15, N2O,Err_N2O, N2O_N15, Err_N2O_N15
#-----------------------------------------------------------------------------------------------------------------------------------------------

obs<-read.csv('Data_1.csv')
names(obs)<-c('time', 'NH4', 'NO3', 'NH4_N15', 'NO3_N15','TON_N15')
obs[,'NH4_N15']<-obs[,'NH4']*obs[,'NH4_N15']/100
obs[,'NO3_N15']<-obs[,'NO3']*obs[,'NO3_N15']/100
obs[,'TON_N15']<-1307*obs[,'TON_N15']/100

times<- obs[,'time']

### MCMC simulation
#----------------------------------------------------------------------------
pp<-pars0*SS;
upper<-lower<-pars0;
lower[]= 0;
upper[]= 1;

system.time({

MCMC<-mdMCMC1(ff=FRAUZcost1,pp=pp,obs=obs,niter=150000,y=y0,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=5000);
MCMC_SS<-mdMCMC(ff=FRAUZcost,pp=pp,obs=obs,niter=150000,y=y0,SS,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=5000);
})

### MCMC output
#------------------------------------------------------------------------------
print(summary(MCMC))
print(summary(MCMC_SS))
### MCMC best parameters
#------------------------------------------------------------------------------
print(MCMC$bestpar)
print(MCMC_SS$bestpar)
###Plot the processes
#-----------------------------------------------------------------------------
win.graph()
cumuplot(as.mcmc(MCMC$pars))

### Plot the state variables
#-----------------------------------------------------------------------------
sim_data<-as.data.frame(FRAUZ(y0,MCMC_SS$bestpar,SS,times));
pv2<-valid_plotMCMC(obs,sim_data);
sim_data<-as.data.frame(FRAUZ1(y0,MCMC$bestpar,times));
pv2<-valid_plotMCMC(obs,sim_data);

### Sensitivity ranges
#------------------------------------------------------------------------------
arg_g=c('NH4','NO3','NH4_N15','NO3_N15','TON_N15','Mib','Hum','NH3','N2O')
sR1 <- sensRange(func = FRAUZ1, parms = MCMC$bestpar, times=times, y=y0, parInput = MCMC$par, dist = "latin",sensvar=arg_g)
	win.graph()
	opar<-par();
	par(mfcol=c(4,5), mar=c(3,3,2,0.5))
	plot(summary(sR1), xlab = "time")
	par(opar);
#-----------------------------------------------------------------------------
dyn.unload('flauz.dll')
