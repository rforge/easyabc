## FUNCTION ABC_rejection: brute-force ABC (Pritchard et al. 1999)
######################################################

ABC_rejection_cluster <-function(model,prior_matrix,nb_simul,seed_count=0,n_cluster=1){
    ## checking errors in the inputs
    if(missing(model)) stop("'model' is missing")
    if(missing(prior_matrix)) stop("'prior_matrix' is missing")
    if(missing(nb_simul)) stop("'nb_simul' is missing")
    if(!is.matrix(prior_matrix) && !is.data.frame(prior_matrix)) stop("'prior_matrix' has to be a matrix or data.frame.")
    if(is.data.frame(prior_matrix)) prior_matrix <- as.matrix(prior_matrix)
    if(dim(prior_matrix)[2]!=2) stop("'prior_matrix' must have two columns.")
    if (nb_simul<1) stop("'nb_simul' must be a number larger than 1.")
    nb_simul=floor(nb_simul)
    if(!is.vector(seed_count)) stop("'seed_count' has to be a number.")
    if(length(seed_count)>1) stop("'seed_count' has to be a number.")
    if (seed_count<0) stop ("'seed_count' has to be a positive number.")
    seed_count=floor(seed_count)
    if(!is.vector(n_cluster)) stop("'n_cluster' has to be a number.")
    if(length(n_cluster)>1) stop("'n_cluster' has to be a number.")
    if (n_cluster<1) stop ("'n_cluster' has to be a positive number.")
    n_cluster=floor(n_cluster)
	library(parallel)

	start = Sys.time()
	
	options(scipen=50)
	cl <- makeCluster(getOption("cl.cores", n_cluster))

	tab_simul_summarystat=NULL
	tab_param=NULL
	list_param=list(NULL)
	npar=floor(nb_simul/n_cluster)
	n_end=nb_simul-(npar*n_cluster)
	for (irun in 1:npar){
	  for (i in 1:n_cluster){
		l=dim(prior_matrix)[1]
		param=array(0,l)
		for (j in 1:l){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		#if (use_seed) { # NB: we force the value use_seed=TRUE
		param=c(n_cluster,(seed_count+i),param) # the first parameter is the number of cores/clusters used
		list_param[[i]]=param
		tab_param=rbind(tab_param,param[3:(l+2)])
	  }
	  seed_count=seed_count+n_cluster
	  list_simul_summarystat=parLapply(cl,list_param,model)
	  for (i in 1:n_cluster){
		tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
	  }
	}
	if (n_end>0){
	  stopCluster(cl)
	  cl <- makeCluster(getOption("cl.cores", 1))
	  list_param=list(NULL)
	  for (i in 1:n_end){
		l=dim(prior_matrix)[1]
		param=array(0,l)
		for (j in 1:l){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		#if (use_seed) { # NB: we force the value use_seed=TRUE
		param=c(n_cluster,(seed_count+i),param) # the first parameter is the number of cores/clusters used
		list_param[[i]]=param
		tab_param=rbind(tab_param,param[3:(l+2)])
	  }
	  seed_count=seed_count+n_end
	  list_simul_summarystat=parLapply(cl,list_param,model)
	  for (i in 1:n_end){
		tab_simul_summarystat=rbind(tab_simul_summarystat,as.numeric(list_simul_summarystat[[i]]))
	  }
    	  stopCluster(cl)
	}
	else{
	  stopCluster(cl)
	}
	options(scipen=0)
	sd_simul=sapply(as.data.frame(tab_simul_summarystat),sd)
list(param=tab_param,stats=tab_simul_summarystat,weights=array(1/nb_simul,nb_simul),stats_normalization=sd_simul,nsim=nb_simul,computime=as.numeric(difftime(Sys.time(), start, unit="secs")))
}
