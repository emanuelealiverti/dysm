
reshape_res = function(out, nn = names(inits), 
		       cl = dimnames(Y)$country,
		       n_years = hyp$Time, 
		       n_countries = hyp$N_countries, 
		       n_clust = hyp$H,  pG = 2) {

	# I primi sono a posto così. Gli altri serve un po' di reshape manuale
	# prendiamo quelli bidimensionali
	nn_red = grep("^drift_|s_", nn,invert=T,value=T)
	res = list()

	# actual Parameters, [nit x countries x years] format
	for(i in 1:length(nn_red)) {
		tmp = out[,grep(paste0("^",nn_red[i],"\\["), colnames(out) )]
		res[[nn_red[i]]] = drop(array(tmp, dim = c(NROW(tmp), ncol(tmp)/n_years,n_years )))
	}

	# multivariate components: drift
	nn_dr = grep("drift", nn,value=T)
	for(i in 1:length(nn_dr)){
		tmp = out[,grep(paste0("^",nn_dr[i],"\\["), colnames(out) )]
		res[[nn_dr[i]]] = drop(array(tmp, dim = c(NROW(tmp), ncol(tmp)/n_countries,n_countries )))
	} 
	nn_s = grep("s_", nn,value=T)
	for(i in 1:length(nn_dr)){
		tmp = out[,grep(paste0("^",nn_s[i],"\\["), colnames(out) )]
		res[[nn_s[i]]] = drop(array(tmp, dim = c(NROW(tmp), ncol(tmp)/n_countries,n_countries )))
	}


	

	return(res)
}

## We use a different parametrisation
add_w = function(pl) {
	w0=w1=w2=array(NA, dim = dim(pl$eta1))
	for(j in 1:dim(pl$eta1)[2]){
	for(k in 1:dim(pl$eta1)[3]){
		tmp = t(apply(cbind(pl$eta1[,j,k], pl$eta2[,j,k]), 1, function(r) eta_to_w(eta1 = r[1], eta2 = r[2])))
		w2[,j,k] = tmp[,3] 
		w1[,j,k] = tmp[,2]
		w0[,j,k] = tmp[,1]
	}
	}
	return(c(pl, list("w0" = w0,"w1" = w1,"w2" = w2)))
}


add_skew_moments = function(pl) {
	# computes sn moments
	a2d = function(a) a / sqrt( 1 + a^2 )
	m = function(xi, omega, d) xi + omega * d * sqrt(2/pi)
	v = function(omega, d) omega^2 * (1 - 2 * d^2 / pi )

	sn_mean=sn_sd=array(NA, dim = dim(pl$xi))
	for(j in 1:dim(pl$xi)[2]){
	for(k in 1:dim(pl$xi)[3]){
		sn_mean[,j,k] = m(pl$xi[,j,k],sqrt(exp(pl$lomega[,j,k])), a2d(pl$alpha[,j,k]))
		sn_sd[,j,k] = sqrt(v(sqrt(exp(pl$lomega[,j,k])), a2d(pl$alpha[,j,k])))
	}}
	return(c(pl, list("sn_mean"  = sn_mean, "sn_sd" = sn_sd)))


}



eta_to_w = function(eta1,eta2) {
	expit = function(x) exp(x) / (1 + exp(x) )
		c(expit(-eta1) * expit(-eta2), expit(eta2) * expit(-eta1), expit(eta1))
}



# This function obtains predictions recursively
get_rec_pred = function(nyears = 10,  pl, hyp)
{
	obj_size = c(NROW(pl$xi),hyp$N_countries, nyears)

	tmp = list()
	tmp$xi     =  array(NA, dim = obj_size)
	tmp$lomega =  array(NA, dim = obj_size)
	tmp$alpha =  array(NA, dim = obj_size)

	tmp$mu     =  array(NA, dim = obj_size)
	tmp$lsigma =  array(NA, dim = obj_size)

	tmp$eta1   =  array(NA, dim = obj_size)
	tmp$eta2   =  array(NA, dim = obj_size)


	for(it in 1:NROW(pl$mu)){
		for(nc in 1:obj_size[2]) { # loop across countries

			tmp$xi[it,nc,1]     = pl$xi[it,nc,hyp$Time] + pl$drift_xi[it,nc]
			tmp$lomega[it,nc,1] = pl$lomega[it,nc,hyp$Time] + pl$drift_lomega[it,nc]
			tmp$alpha[it,nc,1]  = pl$alpha[it,nc,hyp$Time] + pl$drift_alpha[it,nc]

			tmp$mu[it,nc,1]     = pl$mu[it,nc,hyp$Time] + pl$drift_mu[it,nc]
			tmp$lsigma[it,nc,1] = pl$lsigma[it,nc,hyp$Time] + pl$drift_lsigma[it,nc]

			tmp$eta1[it,nc,1]   = pl$eta1[it,nc,hyp$Time] + pl$drift_eta1[it,nc]
			tmp$eta2[it,nc,1]   = pl$eta2[it,nc,hyp$Time] + pl$drift_eta2[it,nc]

			# Now loop in time
			for(tt in 2:obj_size[3]) { # loop in time and get recursive predictions
				tmp$xi[it,nc,tt]     = tmp$xi[it,nc,tt-1] + pl$drift_xi[it,nc]
				tmp$lomega[it,nc,tt] = tmp$lomega[it,nc,tt-1] + pl$drift_lomega[it,nc]
				tmp$alpha[it,nc,tt]  = tmp$alpha[it,nc,tt-1] + pl$drift_alpha[it,nc]

				tmp$mu[it,nc,tt]     = tmp$mu[it,nc,tt-1] + pl$drift_mu[it,nc]
				tmp$lsigma[it,nc,tt] = tmp$lsigma[it,nc,tt-1] + pl$drift_lsigma[it,nc]

				tmp$eta1[it,nc,tt]   = tmp$eta1[it,nc,tt-1] + pl$drift_eta1[it,nc]
				tmp$eta2[it,nc,tt]   = tmp$eta2[it,nc,tt-1] + pl$drift_eta2[it,nc]
			}
			}

		if(it %% 50 == 0) cat(it, " over ", obj_size[1], "\n")
	}
	return(tmp)
}


# This function obtains predictions recursively
get_pred_ppd = function(nyears = 10,  pl, hyp, seed = 1)
{
	set.seed(seed)
	obj_size = c(NROW(pl$xi),hyp$N_countries, nyears)

	tmp = list()
	tmp$xi     =  array(NA, dim = obj_size)
	tmp$lomega =  array(NA, dim = obj_size)
	tmp$alpha =  array(NA, dim = obj_size)

	tmp$mu     =  array(NA, dim = obj_size)
	tmp$lsigma =  array(NA, dim = obj_size)

	tmp$eta1   =  array(NA, dim = obj_size)
	tmp$eta2   =  array(NA, dim = obj_size)

	for(it in 1:NROW(pl$mu)){
		for(nc in 1:obj_size[2]) { # loop across countries
			tmp$xi[it,nc,1]     = rnorm(1, mean = pl$xi[it,nc,hyp$Time] + pl$drift_xi[it,nc], sd         = sqrt(pl$s_xi[it,nc]))
			tmp$lomega[it,nc,1] = rnorm(1, mean = pl$lomega[it,nc,hyp$Time] + pl$drift_lomega[it,nc], sd = sqrt(pl$s_lomega[it,nc]))
			tmp$alpha[it,nc,1]  = rnorm(1, mean = pl$alpha[it,nc,hyp$Time] + pl$drift_alpha[it,nc], sd   = sqrt(pl$s_alpha[it,nc]))

			tmp$mu[it,nc,1]     = rnorm(1, mean = pl$mu[it,nc,hyp$Time] + pl$drift_mu[it,nc], sd         = sqrt(pl$s_mu[it,nc]))
			tmp$lsigma[it,nc,1] = rnorm(1, mean = pl$lsigma[it,nc,hyp$Time] + pl$drift_lsigma[it,nc], sd = sqrt(pl$s_lsigma[it,nc]))

			tmp$eta1[it,nc,1]   = rnorm(1, mean = pl$eta1[it,nc,hyp$Time] + pl$drift_eta1[it,nc], sd     = sqrt(pl$s_eta1[it,nc]))
			tmp$eta2[it,nc,1]   = rnorm(1, mean = pl$eta2[it,nc,hyp$Time] + pl$drift_eta2[it,nc], sd     = sqrt(pl$s_eta2[it,nc]))

			# Now loop in time
			for(tt in 2:obj_size[3]) { # loop in time and get recursive predictions
				tmp$xi[it,nc,tt]     = rnorm(1, mean = tmp$xi[it,nc,tt-1] + pl$drift_xi[it,nc], sd         = sqrt(pl$s_xi[it,nc]))
				tmp$lomega[it,nc,tt] = rnorm(1, mean = tmp$lomega[it,nc,tt-1] + pl$drift_lomega[it,nc], sd = sqrt(pl$s_lomega[it,nc]))
				tmp$alpha[it,nc,tt]  = rnorm(1, mean = tmp$alpha[it,nc,tt-1] + pl$drift_alpha[it,nc], sd   = sqrt(pl$s_alpha[it,nc]))

				tmp$mu[it,nc,tt]     = rnorm(1, mean = tmp$mu[it,nc,tt-1] + pl$drift_mu[it,nc], sd         = sqrt(pl$s_mu[it,nc]))
				tmp$lsigma[it,nc,tt] = rnorm(1, mean = tmp$lsigma[it,nc,tt-1] + pl$drift_lsigma[it,nc], sd = sqrt(pl$s_lsigma[it,nc]))

				tmp$eta1[it,nc,tt]   = rnorm(1, mean = tmp$eta1[it,nc,tt-1] + pl$drift_eta1[it,nc], sd     = sqrt(pl$s_eta1[it,nc]))
				tmp$eta2[it,nc,tt]   = rnorm(1, mean = tmp$eta2[it,nc,tt-1] + pl$drift_eta2[it,nc], sd     = sqrt(pl$s_eta2[it,nc]))
			}

		}
		if(it %% 50 == 0) cat(it, " over ", obj_size[1], "\n")
	}
	return(tmp)
}



get_df_summ = function(pl,pred, countries = dimnames(dd$y)[[1]], 
		       log_scale = T, w = F, 
		       moments = F,
		       func = "mean", ...) {
	df_pl = list()

	for(p in 1:length(countries)){
		cur_c = countries[p]
		tmp = data.frame(x=1:(hyp$Time+dim(pred$xi)[3]))

		if(w) {
			tmp$w0 = c(apply(pl$w0[,p,],2,get(func), ...), apply(pred$w0[,p,],2,get(func), ...))
			tmp$w1 = c(apply(pl$w1[,p,],2,get(func), ...), apply(pred$w1[,p,],2,get(func), ...))
			tmp$w2 = c(apply(pl$w2[,p,],2,get(func), ...), apply(pred$w2[,p,],2,get(func), ...))
		} else {
			tmp$eta1 = c(apply(pl$eta1[,p,],2,get(func), ...), apply(pred$eta1[,p,],2,get(func), ...))
			tmp$eta2 =c( apply(pl$eta2[,p,],2,get(func), ...),  apply(pred$eta2[,p,],2,get(func), ...))
		}

		tmp$mu =c( apply(pl$mu[,p,],2,get(func), ...),  apply(pred$mu[,p,],2,get(func), ...))

		if(log_scale){
			tmp$sigma =c( apply((pl$lsigma[,p,]),2,get(func), ...),  apply((pred$lsigma[,p,]),2,get(func), ...))
			tmp$omega =c( apply((pl$lomega[,p,]),2,get(func), ...),  apply((pred$lomega[,p,]),2,get(func), ...))
		} else {
			tmp$sigma =c( apply(exp(pl$lsigma[,p,]),2,get(func), ...),  apply(exp(pred$lsigma[,p,]),2,get(func), ...))
			tmp$omega =c( apply(exp(pl$lomega[,p,]),2,get(func), ...),  apply(exp(pred$lomega[,p,]),2,get(func), ...))
		}

		tmp$xi =c( apply(pl$xi[,p,],2,get(func), ...),  apply(pred$xi[,p,],2,get(func), ...))
		tmp$alpha =c( apply(pl$alpha[,p,],2,get(func), ...),  apply(pred$alpha[,p,],2,get(func), ...))

		if(moments) {
			tmp$sn_mean =c( apply(pl$sn_mean[,p,],2,get(func), ...),  apply(pred$sn_mean[,p,],2,get(func), ...))
			tmp$sn_sd =c( apply(pl$sn_sd[,p,],2,get(func), ...),  apply(pred$sn_sd[,p,],2,get(func), ...))
			tmp$xi = NULL
			tmp$omega = NULL
		}

		df_pl[[cur_c]] = tmp



	}
	return(df_pl)
}


get_df_ic = function(pl,pred,  
		     countries = dimnames(dd$y)[[1]],
		     log_scale = T, 
		     w = F,
		     moments = F,
		     level = 0.05) {
	df_pl = list()

	for(p in 1:length(countries)){
		cur_c = countries[p]
		tmp = data.frame(x=1:(hyp$Time+dim(pred$xi)[3]))


		if(w) {
			tmp$w0 = c(apply(pl$w0[,p,],2,quantile, level), apply(pred$w0[,p,],2,quantile, level))
			tmp$w1 = c(apply(pl$w1[,p,],2,quantile, level), apply(pred$w1[,p,],2,quantile, level))
			tmp$w2 = c(apply(pl$w2[,p,],2,quantile, level), apply(pred$w2[,p,],2,quantile, level))
		} else {
			tmp$eta1 = c(apply(pl$eta1[,p,],2,quantile, level), apply(pred$eta1[,p,],2,quantile, level))
			tmp$eta2 = c(apply(pl$eta2[,p,],2,quantile, level),  apply(pred$eta2[,p,],2,quantile, level))
		}

		tmp$mu =c( apply(pl$mu[,p,],2,quantile, level),  apply(pred$mu[,p,],2,quantile, level))

		if(log_scale){
			tmp$sigma =c( apply((pl$lsigma[,p,]),2,quantile, level),  apply((pred$lsigma[,p,]),2,quantile, level))
			tmp$omega =c( apply((pl$lomega[,p,]),2,quantile, level),  apply((pred$lomega[,p,]),2,quantile, level))
		} else {
			tmp$sigma =c( apply(exp(pl$lsigma[,p,]),2,quantile, level),  apply(exp(pred$lsigma[,p,]),2,quantile, level))
			tmp$omega =c( apply(exp(pl$lomega[,p,]),2,quantile, level),  apply(exp(pred$lomega[,p,]),2,quantile, level))
		}

		tmp$xi =c( apply(pl$xi[,p,],2,quantile, level),  apply(pred$xi[,p,],2,quantile, level))
		tmp$alpha =c( apply(pl$alpha[,p,],2,quantile, level),  apply(pred$alpha[,p,],2,quantile, level))
		if(moments) {
			tmp$sn_mean =c( apply(pl$sn_mean[,p,],2,quantile, level),  apply(pred$sn_mean[,p,],2,quantile, level))
			tmp$sn_sd =c( apply(pl$sn_sd[,p,],2,quantile, level),  apply(pred$sn_sd[,p,],2,quantile, level))
			tmp$xi = NULL
			tmp$omega = NULL
		}

		df_pl[[cur_c]] = tmp
	}
	return(df_pl)
}




get_df_ic_h = function(pl,pred,  
		       countries = dimnames(dd$y)[[1]],
		       log_scale = T, w = F, moments = F,
		       level = 0.95, upper = T){
	df_pl = list()

	emp_hpd = function(x, level) TeachingDemos::emp.hpd(x, level)[ifelse(upper,2,1)]

	for(p in 1:length(countries)){
		cur_c = countries[p]
		tmp = data.frame(x=1:(hyp$Time+dim(pred$xi)[3]))

		if(w) {
			tmp$w0 = c(apply(pl$w0[,p,],2,emp_hpd, level), apply(pred$w0[,p,],2,emp_hpd, level))
			tmp$w1 = c(apply(pl$w1[,p,],2,emp_hpd, level), apply(pred$w1[,p,],2,emp_hpd, level))
			tmp$w2 = c(apply(pl$w2[,p,],2,emp_hpd, level), apply(pred$w2[,p,],2,emp_hpd, level))
		} else {
			tmp$eta1 = c(apply(pl$eta1[,p,],2,emp_hpd, level), apply(pred$eta1[,p,],2,emp_hpd, level))
			tmp$eta2 = c(apply(pl$eta2[,p,],2,emp_hpd, level),  apply(pred$eta2[,p,],2,emp_hpd, level))
		}


		tmp$mu =c( apply(pl$mu[,p,],2,emp_hpd, level),  apply(pred$mu[,p,],2,emp_hpd, level))

		if(log_scale){
			tmp$sigma =c( apply((pl$lsigma[,p,]),2,emp_hpd, level),  apply((pred$lsigma[,p,]),2,emp_hpd, level))
			tmp$omega =c( apply((pl$lomega[,p,]),2,emp_hpd, level),  apply((pred$lomega[,p,]),2,emp_hpd, level))
		} else {
			tmp$sigma =c( apply(exp(pl$lsigma[,p,]),2,emp_hpd, level),  apply(exp(pred$lsigma[,p,]),2,emp_hpd, level))
			tmp$omega =c( apply(exp(pl$lomega[,p,]),2,emp_hpd, level),  apply(exp(pred$lomega[,p,]),2,emp_hpd, level))
		}

		tmp$xi =c( apply(pl$xi[,p,],2,emp_hpd, level),  apply(pred$xi[,p,],2,emp_hpd, level))
		tmp$alpha =c( apply(pl$alpha[,p,],2,emp_hpd, level),  apply(pred$alpha[,p,],2,emp_hpd, level))

		if(moments) {
			tmp$sn_mean =c( apply(pl$sn_mean[,p,],2,emp_hpd, level),  apply(pred$sn_mean[,p,],2,emp_hpd, level))
			tmp$sn_sd =c( apply(pl$sn_sd[,p,],2,emp_hpd, level),  apply(pred$sn_sd[,p,],2,emp_hpd, level))
			tmp$xi = NULL
			tmp$omega = NULL
		}

		df_pl[[cur_c]] = tmp
	}
	return(df_pl)
}




