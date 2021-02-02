library(nimble)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Link function corresponding to a discretized mixture density function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(!file.exists("skew_cdf.o")){
if(!fPIC){
	system("g++ -c -Wall -I ./ skew_cdf.cpp")
} else {
	system("g++ -c -Wall -fPIC -I ./ skew_cdf.cpp") # under some architectures
}
}


diffc <- nimbleExternalCall( function(out = double(1),  # output
				     x = double(1), # input
				n = integer(0)
				){},
			    Cfun =  'diffc',
		    headerFile = file.path(getwd(),"./skew_cdf.hpp"),
			    returnType = void(), 
			    oFile = file.path(getwd(),"./skew_cdf.o"))

diffc_vec = nimbleFunction(run = function(x = double(1)) {
				   val <- numeric(length(x)) 
				   diffc(val,x,length(x))
				   return(val)
				            returnType(double(1))

})

psnc <- nimbleExternalCall( function(out = double(1),  # output
				     x = double(1), # input
				      xi = double(0),
				      omega = double(0),
				      alpha = double(0),
				n = integer(0)
				){},
			    Cfun =  'psnc',
		    headerFile = file.path(getwd(),"./skew_cdf.hpp"),
			    returnType = void(), 
			    oFile = file.path(getwd(),"./skew_cdf.o"))

psnc_vec = nimbleFunction(run = function(x = double(1), 
					 xi = double(0), 
					 omega = double(0), 
					 alpha = double(0) ) {
				   val <- numeric(length(x)) 
				   psnc(val,x, xi,omega,alpha,length(x))
				   return(val)
				            returnType(double(1))
})



link_func = nimbleFunction(
			   run = function(x=double(1), 
					  w = double(1,3), 

					  mu = double(0),
					  lsigma = double(0),

					  xi = double(0),
					  lomega = double(0),
					  alpha = double(0)
					  ){
				   returnType(double(1))

				   # Skek normal cdf (check for positive contraints)
				   tmp3 <- psnc_vec(x, xi,exp(lomega),alpha) 
				   tmp2 <- pnorm(x, mu, sd = exp(lsigma))
				   tmp <- w[2] * diffc_vec(tmp2) + w[3] * diffc_vec(tmp3)
				   tmp[1] <- w[1]


				   #return(tmp)
				   return(tmp/sum(tmp))
			   }
)
