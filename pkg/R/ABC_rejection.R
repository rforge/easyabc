## FUNCTION ABC_rejection: brute-force ABC (Pritchard et al. 1999)
######################################################
ABC_rejection<-function(model,prior,nb_simul,use_seed=FALSE,seed_count=0,n_cluster=1,verbose=FALSE,progress_bar=FALSE,summary_stat_target,tol){
    ## checking errors in the inputs
    if(missing(model)) stop("'model' is missing")
    if(missing(prior)) stop("'prior' is missing")
    if(!is.list(prior)) stop("'prior' has to be a list")
    l=length(prior)
    for (i in 1:l){
    	if(!any(prior[[i]][1] == c("unif", "normal", "lognormal", "exponential"))) {
        	stop("Prior distribution type must be unif, normal, lognormal or exponential")
    	}
	if (prior[[i]][1]=="exponential"){
		if (length(prior[[i]])<2){
			stop(paste("Incomplete prior information for parameter ",i,sep=""))
		}
	}
	else{
		if (length(prior[[i]])<3){
			stop(paste("Incomplete prior information for parameter ",i,sep=""))
		}
	}
    }
    if(missing(nb_simul)) stop("'nb_simul' is missing")
    if (nb_simul<1) stop("'nb_simul' must be a number larger than 1")
    if(!is.logical(use_seed)) stop("'use_seed' has to be boolean")
    if(!is.vector(seed_count)) stop("'seed_count' has to be a number")
    if(length(seed_count)>1) stop("'seed_count' has to be a number")
    if (seed_count<0) stop ("'seed_count' has to be a positive number")
    if(!is.vector(n_cluster)) stop("'n_cluster' has to be a number.")
    if(length(n_cluster)>1) stop("'n_cluster' has to be a number.")
    if (n_cluster<1) stop ("'n_cluster' has to be a positive number.")
    n_cluster=floor(n_cluster)
    if(!is.logical(verbose)) stop("'verbose' has to be boolean")
    if(!is.logical(progress_bar)) stop("'progress_bar' has to be boolean")

    nb_simul=floor(nb_simul)
    seed_count=floor(seed_count)
    rejection=NULL
    if (n_cluster==1){
    	rejection=.ABC_rejection(model,prior_matrix,nb_simul,use_seed,seed_count,verbose,progress_bar)
    }
    else{
	if (use_seed==FALSE){
		stop("For parallel implementations, you must specify the option 'use_seed=TRUE' and modify your model accordingly - see the package's vignette for more details.")
	}
	rejection=.ABC_rejection_cluster(model,prior_matrix,nb_simul,seed_count,n_cluster,verbose)
    }

    res=NULL
    if (missing(summary_stat_target)){
	res=list(param=rejection$param, stats=rejection$stats, weights=rejection$weights, stats_normalization=rejection$stats_normalization, nsim=rejection$nsim, computime=rejection$computime)
    }
    else{
	if (missing(tol)){
		stop("'tol' is missing")
	}
	rej=abc(summary_stat_target,rejection$param,rejection$stats,tol,method="rejection")
	nr=dim(rej$unadj.values)[1]
	res=list(param=rej$unadj.values, stats=rej$ss, weights=array(1/nr,nr), stats_normalization=rejection$stats_normalization, nsim=rejection$nsim, nrec=nr, computime=rejection$computime)
    }
res
}


