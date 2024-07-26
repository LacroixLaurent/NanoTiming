### rescaling function
myscaling0 <- function(x,infq=0.005,supq=0.995,...)
	{
		upper <- quantile(x,supq)
		lower <- quantile(x,infq)
		output <- (x-lower)/(upper-lower)
		return(output)
	}
	
# this rescaling function put the signal to 0 for the quantile corresponding to infq and to 1 for the signal corresponding to supq.