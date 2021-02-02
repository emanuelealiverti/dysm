source("./nimble_utilities.r")
mod_multi_country_mixture <- nimbleCode({ 


	for(p in 1:N_countries){


		s_xi[p]     ~ dinvgamma(1,1)
		s_lomega[p] ~ dinvgamma(1,1)
		s_alpha[p]  ~ dinvgamma(1,1)
		s_mu[p]     ~ dinvgamma(1,1)
		s_lsigma[p] ~ dinvgamma(1,1)
		s_eta1[p]   ~ dinvgamma(1,1)
		s_eta2[p]   ~ dinvgamma(1,1)



		# Drift parameters
		drift_xi[p] ~ dnorm(mean = 0, sd = 1)
		drift_lomega[p] ~ dnorm(mean = 0, sd = 1)
		drift_alpha[p] ~ dnorm(mean = 0, sd = 1)

		drift_mu[p] ~ dnorm(mean = 0, sd = 1)
		drift_lsigma[p] ~ dnorm(mean = 0, sd = 1)

		drift_eta1[p] ~ dnorm(mean = 0, sd = 1)
		drift_eta2[p] ~ dnorm(mean = 0, sd = 1)


		m_xi_0[p] ~ dnorm( mean = m_xi, sd = 1)
		m_lomega_0[p] ~ dnorm( mean = m_lomega, sd = 1)
		m_alpha_0[p] ~ dnorm( mean = m_alpha, sd = 1)
		m_mu_0[p] ~ dnorm( mean = m_mu, sd = 1)
		m_lsigma_0[p] ~ dnorm( mean = m_lsigma, sd = 1)
		m_eta1_0[p] ~ dnorm( mean = m_eta1, sd = 1)
		m_eta2_0[p] ~ dnorm( mean = m_eta2, sd = 1)

		# random walk. First observation
		xi[p,1] ~ dnorm(mean     =  m_xi_0[p], var = s_xi[p])
		lomega[p,1] ~ dnorm(mean =  m_lomega_0[p], var = s_lomega[p])
		alpha[p,1] ~ dnorm(mean  =  m_alpha_0[p], var = s_alpha[p])

		mu[p,1] ~ dnorm(mean     =  m_mu_0[p], var = s_mu[p])
		lsigma[p,1] ~ dnorm(mean =  m_lsigma_0[p], var = s_lsigma[p])

		eta1[p,1] ~ dnorm(mean = m_eta1_0[p], var = s_eta1[p])
		eta2[p,1] ~ dnorm(mean = m_eta2_0[p], var = s_eta2[p])
	}



	#++++++++++++++++++++++++++++++++++++++++++++++++++
	# Time recursion: multivariate random walk + drift
	#++++++++++++++++++++++++++++++++++++++++++++++++++

	for(p in 1:N_countries){
		for (tt in 2:Time){
			xi[p,tt] ~ dnorm(mean = drift_xi[p] + xi[p, tt-1], var = s_xi[p])
			lomega[p,tt] ~ dnorm(mean = drift_lomega[p] + lomega[p, tt-1], var = s_lomega[p])
			alpha[p,tt] ~ dnorm(mean = drift_alpha[p] + alpha[p, tt-1], var = s_alpha[p])

			mu[p,tt] ~ dnorm(mean = drift_mu[p] + mu[p, tt-1], var = s_mu[p])
			lsigma[p,tt] ~ dnorm(mean = drift_lsigma[p] + lsigma[p, tt-1], var = s_lsigma[p])

			eta1[p,tt] ~ dnorm(mean = drift_eta1[p] + eta1[p, tt-1], var = s_eta1[p])
			eta2[p,tt] ~ dnorm(mean = drift_eta2[p] + eta2[p, tt-1], var = s_eta2[p])

		}

	}



	# multinomial likelihood
	for(p in 1:N_countries){

		for(tt in 1:Time) {
			# We use a different parametrization here, cumulative logit
			# works better in practice
			w[p,3,tt] <- expit(eta1[p,tt] )
			w[p,2,tt] <- expit(eta2[p,tt]) * expit(-eta1[p,tt])
			w[p,1,tt] <- 1 - w[p,2,tt] - w[p,3,tt]
			#w[p,1,tt] <- expit(-eta2[p,tt]) * expit(-eta1[p,tt])

			probs[p,1:n,tt] <- link_func(x[1:n],w[p,1:3,tt],
						     mu[p,tt],lsigma[p,tt],
						     xi[p,tt],lomega[p,tt],alpha[p,tt])
			y[p,1:n,tt] ~ dmulti(probs[p,1:n,tt],N[p,tt])

		}


	}
})
