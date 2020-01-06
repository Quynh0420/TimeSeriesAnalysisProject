compare.ICs.GARCH<-function(models_list)
{ 
	n_ <- length(models_list)
	
	for(i in 1:n_)
		{
		ICs_ <- data.frame(t(get(models_list[i])@fit$ics))
		ICs_$model <- models_list[i]
		if(i==1) ICs <- ICs_ else ICs <- rbind(ICs, ICs_)
		}
	
    mins_ <- sapply(ICs[,1:(ncol(ICs)-1)],
                function(x) which(x == min(x)))
	
	return(list(ICs = ICs, which.min = mins_))
}