ABC_rejection <-
function(model,prior_matrix,nb_simul,use_seed=TRUE,seed_count=0){
	tab_simul_summarystat=NULL
	tab_param=NULL
  
	start = Sys.time()
	
	
  # progress bar
	pb <- .progressBar(width=50)
  for (i in 1:nb_simul){
		l=dim(prior_matrix)[1]
		param=array(0,l)
		for (j in 1:l){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		if (use_seed) {
			param=c((seed_count+i),param)
		}
		simul_summarystat=model(param)
		tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
		if (use_seed) {
			tab_param=rbind(tab_param,param[2:(l+1)])
		}
		else{
			tab_param=rbind(tab_param,param)
		}
    
    # for progressbar message and time evaluation 
		duration = duration + difftime(Sys.time(), start, unit="secs")
		text="";
		if (i==length(simus)) {
		   text = paste("Completed  in",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
		} else {
		 text = paste("− Time elapsed:",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"− Estimated time remaining:",format(.POSIXct(duration/i*(length(simus)-i), tz="GMT"), "%H:%M:%S"));
		}
		.updateProgressBar(pb, i/length(simus), text)
  }
	close(pb)
	cbind(tab_param,tab_simul_summarystat)
}

