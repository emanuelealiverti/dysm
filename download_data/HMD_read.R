                                        #hmd functions to get life table mx and dx and actual dx, functions taken from package demography
hmd.mx <- function (country, username, password, label = country)
    {
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", 
        "Mx_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    mx <- try(utils::read.table(con, skip = 2, header = TRUE, 
        na.strings = "."), TRUE)
    close(con)
    if (class(mx) == "try-error") 
        stop("Connection error at www.mortality.org. Please check username, password and country label.")
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", 
        "Exposures_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    pop <- try(utils::read.table(con, skip = 2, header = TRUE, 
        na.strings = "."), TRUE)
    close(con)
    if (class(pop) == "try-error") 
        stop("Exposures file not found at www.mortality.org")
    obj <- list(type = "mortality", label = label, lambda = 0)
    obj$year <- sort(unique(mx[, 1]))
    n <- length(obj$year)
    m <- length(unique(mx[, 2]))
    obj$age <- mx[1:m, 2]
    mnames <- names(mx)[-c(1, 2)]
    n.mort <- length(mnames)
    obj$rate <- obj$pop <- list()
    for (i in 1:n.mort) {
        obj$rate[[i]] <- matrix(mx[, i + 2], nrow = m, ncol = n)
        obj$rate[[i]][obj$rate[[i]] < 0] <- NA
        obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
        obj$pop[[i]][obj$pop[[i]] < 0] <- NA
        dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, 
            obj$year)
    }
    names(obj$pop) = names(obj$rate) <- tolower(mnames)
    obj$age <- as.numeric(as.character(obj$age))
    if (is.na(obj$age[m])) 
        obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
    return(structure(obj, class = "demogdata"))
    }

#Reads Number of deaths from HMD
hmd.deaths <- function (country, username, password, label = country)
    {
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", 
        "Deaths_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    deaths <- try(utils::read.table(con, skip = 2, header = TRUE, 
        na.strings = "."), TRUE)
    close(con)
    if (class(deaths) == "try-error") 
        stop("Connection error at www.mortality.org. Please check username, password and country label.")
     return(deaths)
    }
                                        #Transform mortality rates into death counts of life table
mx2dx <- function(mx,age,sex="F"){
    ageint <- age[-1]-age[-length(age)]
    ax <- rep(NA,length(age))
    ax[1:length(age)-1] <- ageint/2
    ax[1] <- as.numeric(sex=="F")*as.numeric(mx[1]>=0.107)*0.350+as.numeric(sex=="M")*as.numeric(mx[1]>=0.107)*0.330+
        as.numeric(sex=="F")*as.numeric(mx[1]<0.107)*(0.053+2.8*mx[1])+as.numeric(sex=="M")*as.numeric(mx[1]<0.107)*(0.045+2.684*mx[1])
    ax[length(age)] <- 1/mx[length(age)]
    qx <- rep(NA,length(age))
    qx[1:length(age)-1] <- ageint*mx[1:length(age)-1]/(1+(ageint-ax[1:length(age)-1])*mx[1:length(age)-1])
    qx[length(age)] <- 1
    lx <- rep(NA,length(age))
    lx <- round(c(100000,100000*cumprod(1-qx[-length(age)])))
    dx <- rep(NA,length(age))
    dx[1:length(age)-1] <- lx[-length(age)]-lx[-1]
    dx[length(age)] <- lx[length(age)]
    dx[is.na(dx)] <- 0
    return(dx)
}
