## FUNCTION ABC_rejection: brute-force ABC (Pritchard et al. 1999)
######################################################

ABC_rejection<-function(model,prior_matrix,nb_simul,use_seed=TRUE,seed_count=0){

    ## checking errors in the inputs
    if(missing(model)) stop("'model' is missing")
    if(missing(prior_matrix)) stop("'prior_matrix' is missing")
    if(missing(nb_simul)) stop("'nb_simul' is missing")
    if(!is.matrix(prior_matrix) && !is.data.frame(prior_matrix)) 			  stop("'prior_matrix' has to be a matrix or data.frame")
    if(is.data.frame(prior_matrix)) prior_matrix <- as.matrix(prior_matrix)
    if(dim(prior_matrix)[2]!=2) stop("'prior_matrix' must have two columns")
    if (nb_simul<1) stop("'nb_simul' must be a number larger than 1")
    if(!is.logical(use_seed)) stop("'use_seed' has to be boolean")
    if(!is.vector(seed_count)) stop("'seed_count' has to be a number")
    if(length(seed_count)>1) stop("'seed_count' has to be a number")
    if (seed_count<0) stop ("'seed_count' has to be a positive number")

    nb_simul=floor(nb_simul)
    seed_count=floor(seed_count)
    options(scipen=50)
    tab_simul_summarystat=NULL
    tab_param=NULL
    
    start = Sys.time()
    # progress bar
    pb <- .progressBar(width=50)
    duration = 0;
    for (i in 1:nb_simul) {
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
	else {
		tab_param=rbind(tab_param,param)
	}

	# for progressbar message and time evaluation
	duration = difftime(Sys.time(), start, unit="secs")
	text = "";
	if (i == nb_simul) {
	    text = paste("Completed  in", format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
	} 
	else {
	  text = paste("Time elapsed:",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"Estimated time remaining:",format(.POSIXct(duration/i*(nb_simul-i), tz="GMT"), "%H:%M:%S"));
	}
	.updateProgressBar(pb, i/nb_simul, text)
    }
    close(pb)
    options(scipen=0)
    sd_simul=sapply(as.data.frame(tab_simul_summarystat), sd)
    
list(param=tab_param, stats=tab_simul_summarystat, weights=array(1/nb_simul,nb_simul), stats_normalization=sd_simul, nsim=nb_simul, computime=as.numeric(difftime(Sys.time(), start, unit="secs")))
}

