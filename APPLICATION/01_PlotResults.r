# This script loads the results and reproudce fig 5 of the paper

rm(list=ls())
dir()
Npred = 10
p_fold = "./plots/"
dir.create(p_fold,showWarnings = F)

# change here if you want only one
sexes = c("MALE", "FEMALE")
trace_plots = T

pl_list = list()
for (wh in sexes){
	#wh = "FEMALE"



	fPIC =F
	source("./nimble_model_integral.r")
	source("./post_proc_util.r")

	load(sprintf("./mcmc_res%s.RData",wh))
	#pred_c = pred

	#pred = get_rec_pred(10,pars,hyp)
	Y = pars$data
	hyp = list(N_countries = dim(Y)[1], 
		   Time = dim(Y)[3], 
		   n = 111, x = seq(0,110)+0.5)

	# Evaluate difference
	pars = add_w(pars)
	(cc = dimnames(pars$data)[[1]])

	pred = get_pred_ppd(Npred,pars,hyp)
	cc = dimnames(Y)[[1]]
	years = dimnames(Y)[[3]]

	require(ggplot2)
	require(latex2exp)
	require(plyr)

	parsw = add_w(pars)
	predw = add_w(pred)

	#+++++++++++++++++++++
	# Focus now on moments
	#+++++++++++++++++++++
	parsw = add_skew_moments(parsw)
	predw = add_skew_moments(predw)

	df_pl = get_df_summ(parsw,predw, countries = cc, func = "median",  w = T,log=F,moments = T )
	df_qq = get_df_ic_h(parsw,predw, countries = cc, level = .8, w = T, log = F, upper = F, moments = T)
	df_QQ = get_df_ic_h(parsw,predw, countries = cc, level = .8, w = T, log = F, moments = T)
	{
		df_comb = reshape2::melt(df_pl,id.var = c("x"))
		df_comb_QQ = reshape2::melt(df_QQ,id.var = c("x"))
		df_comb_QQ$q = "Q"
		df_comb_qq = reshape2::melt(df_qq,id.var = c("x"))
		df_comb_qq$q = "q"

		from_strings = c("w0", "w1",    "w2", "mu" ,  "sigma" , "sn_mean", "sn_sd", "alpha")
		to_strings = c(TeX("$\\pi_{0t}$"), TeX("$\\pi_{1t}$"), TeX("$\\pi_{2t}$"), 
			       TeX("$\\mu_{t}$"), TeX("$\\sigma_{t}$"), 
			       TeX("$\\tilde{\\mu}^{SN}_{t}$"), 
			       TeX("$\\tilde{\\sigma}^{SN}_{t}$"), 
			       TeX("$\\alpha_{t}$"))

		df_comb$variable = mapvalues(df_comb$variable, 
					     from = from_strings,
					     to = to_strings)
		df_comb$L1 = mapvalues(df_comb$L1, from = c("DEUTW", "FRATNP", "GBRTENW"), to = c("DEU", "FRA", "GBR"))


		df_comb_qq$variable = mapvalues(df_comb_qq$variable, 
						from = from_strings,
						to = to_strings)
		df_comb_qq$L1 = mapvalues(df_comb_qq$L1,from = c("DEUTW", "FRATNP", "GBRTENW"), to = c("DEU", "FRA", "GBR"))


		df_comb_QQ$variable = mapvalues(df_comb_QQ$variable, 
						from = from_strings,
						to = to_strings)
		df_comb_QQ$L1 = mapvalues(df_comb_QQ$L1, from = c("DEUTW", "FRATNP", "GBRTENW"), to = c("DEU", "FRA", "GBR"))


		require(dplyr)
		data_test = df_comb %>% 
			group_by(variable) %>% summarise(m = min(value), M = max(value))

		tmp = years
		(years = c(years[-length(years)], seq(as.numeric(years[length(years)]), len = dim(pred$xi)[3]+1, by =  diff(as.numeric(years))[1])))
		df_comb$x = as.numeric(years[df_comb$x])
		df_comb_QQ$x = as.numeric(years[df_comb_QQ$x])
		df_comb_qq$x = as.numeric(years[df_comb_qq$x])

		df_comb$sex = wh
		df_comb_QQ$sex = wh
		df_comb_qq$sex = wh



		smooth_p = F # if you want to smooth results
		{
			plw = 
				ggplot(df_comb) + 
				facet_wrap(~sex+variable,scales="free_y",nrow=2,labeller = label_parsed) + 
				geom_vline(xintercept = as.numeric(years[hyp$Time]),lty = 3)+

				{if(smooth_p) list(
						   geom_smooth(data=df_comb_QQ,  aes(x,y=value,group=L1, color = L1),lty=3,show.legend=F,se=F) ,
						   geom_smooth(data=df_comb_qq,  aes(x,y=value,group=L1, color = L1),lty=3,show.legend=F,se=F),
						   geom_smooth(aes(x,y=value,group=L1,color=L1),se = F)
				)
						} +
						{if(!smooth_p) list(
								    geom_line(data=df_comb_QQ,  aes(x,y=value,group=L1, color = L1),lty=3,show.legend=F) ,
								    geom_line(data=df_comb_qq,  aes(x,y=value,group=L1, color = L1),lty=3,show.legend=F),
								    geom_line(aes(x,y=value,group=L1,color=L1))
						)
				} +
					#scale_color_manual(values=cbbPalette)+
					theme_bw(base_size = 16) +
					xlab("") + ylab("")+
					scale_x_continuous(breaks = seq(as.numeric(years)[1], as.numeric(tail(years,1)), by = 5)) +
					#scale_x_continuous(breaks = as.numeric(years)[round(seq(1,length(years),by = 5))] ) +
					theme(axis.text.x = element_text(angle = 90,vjust = 1), 
					      strip.text.x = element_text(size=16), 
					      legend.position = ifelse(wh == "MALE", "none", "bottom"),
					      legend.text = element_text(size=12),
					      legend.box = "horizontal",
					      legend.background = element_rect(fill="lightgray", size=0.5, linetype="solid"),
					      legend.title = element_blank(),
					      legend.key = element_rect(colour = 'black')) + 
guides(col = guide_legend(nrow = 2))
		}
	}
	pl_list = c(pl_list, list(plw))
}

library(gtable)
library(ggplot2)

combo_plot = rbind(ggplotGrob(pl_list[[1]] + theme(strip.text.x = element_text(size=16) ) ), 
		   ggplotGrob(pl_list[[2]] + theme(strip.text.x = element_text(size=16) )))


plot(combo_plot)

ggsave(combo_plot, file = "forecast_comb.png", width = 20, height = 17)

