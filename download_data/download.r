# 1] necessario creare un account su https://www.mortality.org/mp/auth.pl

user = "aa@bb.com"
psw = "top_secret_password"

source("./HMD_read.R")

# Tutti i labels sono disponibili https://www.mortality.org/hmd/
# Ad esempio, sud europa + nonrd europa
paesi = c("ITA", "ESP", "PRT", "DNK", "FIN", "SWE", "NOR")

## Scaricare numero totale di decessi e tassi, Può metterci alcuni secondo, a seconda del numero di paesi.
# Se sei sicuro dei label etc vai diretto, sennò meglio con try

# all_data = lapply(paesi, function(pp) list("deaths" =  hmd.deaths(country = pp, username = user,password = psw ), "mx"     =  hmd.mx(country = pp, username = user ,password = psw)))

get_dd = function(pp) list("deaths" =  hmd.deaths(country = pp, username = user,password = psw ), "mx"     =  hmd.mx(country = pp, username = user ,password = psw))

all_data = list()
for(pp in 1:length(paesi)) {
	all_data[[pp]] = c(all_data, try(get_dd(paesi[pp])))
	cat(paesi[pp], " finito \n")
}
names(all_data) = paesi

## cosa ci interessa?
mx$year #anni disponibili
mx$age #classi di età [0 - 110]

# In generale gli anni disponibili non sono consistenti
lapply(all_data, function(pp) range(pp$mx$year))

# Posso prenderli qui dentro, in realtà ha senso anche lavorare su intervalli più piccoli
(interval = c(max(unlist(lapply(all_data, function(pp) min(pp$mx$year)))), min(unlist(lapply(all_data, function(pp) max(pp$mx$year))))))

dx_all = array(NA,dim = c(length(paesi), 1 + interval[2]-interval[1], 2, 111), 
	       dimnames = list("Paesi" = paesi, 
			       "years" = as.character(interval[1]:interval[2]), 
			       "sex" = c("F", "M"), 
			       "age" = as.character(0:110)))

for( pp in 1:length(paesi) ){
	for( curr_year in 1:length(year_grid) ) {
		#pp = 1
		#curr_year = 3
		# Prendi l'anno giusto, diverso per ogni paese
		y = which(year_grid[curr_year] == all_data[[pp]]$mx$year)
		dx_all[pp,curr_year,1,] = mx2dx( all_data[[pp]]$mx$rate$female[,y], age=0:110)/1e5
		dx_all[pp,curr_year,2,] = mx2dx( all_data[[pp]]$mx$rate$male[,y],age=0:110 )/1e5
	}
}
# La curva per il primo paese, primo anno considerato, donne,
plot(dx_all[1,1,1,])
# La curva per il primo paese, primo anno considerato, uomini
plot(dx_all[3,1,2,])
save(dx_all, file = "DX.RData")
