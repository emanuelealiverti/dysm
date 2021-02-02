# This script illustrate the use of DYSM in modelling and forecasting developed countries mortality.
# In order to reproduce results in the paper, you need to download data from HMD creating an account. See the  download_data folder for further instructions

# This script focuses on FEMALE population.
# The Rdata files mcmc_resFEMALE.RData and mcmc_resMALE.RData store the mcmc output
rm(list=ls())

data_fol = "../download_data/"

source(paste0(data_fol, "HMD_read.R"))
cc=sort(c("AUS","DEUTW", "CHE","GBRTENW",    "DNK",    "FIN",    "FRATNP", "ITA" ,   "NLD", "NOR",  "ESP",  "SWE" ))
cc

# LOAD data. Check the content of ../download_data for instructions
ddl = list()
for(i in 1:length(cc)){
	load(file= sprintf(paste0(data_fol, "%s.RData"),cc[i]))
	ddl[[i]] = list("deaths" = deaths,"mx" = mx$rate$female) #change here for Male pop
}
names(ddl) = cc
(Ncount = length(ddl))

Ny = 20 # number of time points
(years= as.character(seq(from  = 2017 - 2*(Ny-1), to = 2017, length.out = Ny))) # You can model year by year, even years, or longer periods

# Combine dx with D
Y = array(NA, c(Ncount,111,Ny),dimnames = list("country" = cc,"ages" = as.character(0:110), "years" = years))
Dtot = array(NA, c(Ncount,Ny),dimnames = list("country" = cc, "years" = as.character(1:Ny)))
for (p in 1:Ncount){
	for(tt in 1:length(years)) {
		(wy = which(years[tt] == (colnames(ddl[[p]]$mx))))
		# calcolo dx
		Dtot[p,tt] = sum(ddl[[p]]$deaths$Female[(ddl[[p]]$deaths$Year == years[tt])])  #Change here for male 
		Y[p,,tt] = mx2dx(ddl[[p]]$mx[,wy], age = 0:110)
		Y[p,,tt] = pmax( round( (Y[p,,tt] / sum (Y[p,,tt])) * Dtot[p,tt]), 0)
	}
}

dd = list(y = Y)
(dd$N = apply(dd$y,c(1,3),sum))
hyp = list(N_countries = dim(dd$y)[1], 
	   Time = dim(dd$y)[3], 
	   n = 111, x = seq(0,110)+0.5)


# Informative priors
hyp$m_mu = 30
hyp$m_xi = 80
hyp$m_alpha = -2
hyp$m_lsigma = 1
hyp$m_lomega = 3
hyp$m_eta1 = 6
hyp$m_eta2 = 1

#++++++++++++++++++++++++++++++++++
# Smart initialisation really helps
#++++++++++++++++++++++++++++++++++
inits = list(mu = rep(hyp$m_mu,Ny), 
	     lsigma = rep(hyp$m_lsigma,Ny), 
	     xi = rep(hyp$m_xi,Ny), 
	     lomega = rep(hyp$m_lomega,Ny),
	     alpha = rep(hyp$m_alpha,Ny),
	     eta1 = rep(hyp$m_eta1,Ny),
	     eta2 = rep(hyp$m_eta2,Ny)
)
inits = lapply(inits, function(x) matrix(x,nrow= hyp$N_countries, ncol = hyp$Time))

inits$drift_mu     =  rep(0, hyp$N_countries)
inits$drift_lsigma =  rep(0, hyp$N_countries)
inits$drift_xi     =  rep(0, hyp$N_countries)
inits$drift_lomega =  rep(0, hyp$N_countries)
inits$drift_alpha  =  rep(0, hyp$N_countries)
inits$drift_eta1   =  rep(0, hyp$N_countries)
inits$drift_eta2   =  rep(0, hyp$N_countries)

inits$s_mu     =  rep(2, hyp$N_countries)
inits$s_lsigma =  rep(2, hyp$N_countries)
inits$s_xi     =  rep(2, hyp$N_countries)
inits$s_lomega =  rep(2, hyp$N_countries)
inits$s_alpha  =  rep(2, hyp$N_countries)
inits$s_eta1   =  rep(2, hyp$N_countries)
inits$s_eta2   =  rep(2, hyp$N_countries)


inits$m_mu_0     =  rep(inits$m_mu, hyp$N_countries)
inits$m_xi_0     =  rep(inits$m_xi, hyp$N_countries)
inits$m_alpha_0  =  rep(inits$m_alpha, hyp$N_countries)
inits$m_lsigma_0 = rep(inits$m_lsigma, hyp$N_countries)
inits$m_lomega_0 = rep(inits$m_lomega, hyp$N_countries)
inits$m_eta1_0   = rep(inits$m_eta1, hyp$N_countries)
inits$m_eta2_0   = rep(inits$m_eta2, hyp$N_countries)




# This depends on your architecture
fPIC = F
source("./nimble_model_integral.r")
source("./post_proc_util.r")
nimbleOptions(stop_after_processing_model_code = F)
dynModel = nimbleModel(code = mod_multi_country_mixture, 
		       constants = hyp, 
		       data = dd, 
		       inits = inits)

cat("============\n[compiling]\n============\n")
CdynModel = compileNimble(dynModel,  showCompilerOutput = TRUE)
cat("============\n[building]\n============\n")
dynMCMC = buildMCMC(dynModel,  monitors = c(names(inits)), print=T)
cat("=========\n[merging]\n============\n")
CdynMCMC = compileNimble(dynMCMC, project = dynModel, showCompilerOutput = TRUE)
CdynModel$initializeInfo()

set.seed(42)
mcmc.out = runMCMC(CdynMCMC, niter = 205000,nburnin = 5000,thin = 50)
pars = reshape_res(mcmc.out,nn = c(names(inits)))
pars = reshape_res(mcmc.out,nn = c(names(inits)))

