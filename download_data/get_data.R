# NB you need to replace user and pw with your credentials
# we provide, for convenience, a list of available countries
source("./HMD_read.R")
cc = read.csv("./countries.csv")
cc$Directory[cc$First_year > 1960]
cc = cc[cc$First_year < 1980,]
(count = sort(unique(cc$Directory)))

## 
get_and_save = function(what = "ITA"){
	deaths = hmd.deaths(country = what, username = "user",password = "pw")
	mx = hmd.mx(country = what, username = "user",password = "pw")
	save(list = c("mx", "deaths"), file = paste0(what,".RData"))
}
count = c("AUT", "DEU",    "CHE","GBRTENW",    "DNK",    "FIN",    "FRATNP", "ITA" ,   "NLD", "NOR",  "ESP",  "SWE" )
sapply(count, function(p) try(get_and_save(what = p)))

# if R < 3.5.0

#for(what in count) {
#load(file = paste0(what,".RData"))
	#save(list = c("mx", "deaths"), file = paste0(what,"v2.RData"), version = 2)
#}
