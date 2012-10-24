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

    rejection = .ABC_rejection_internal(model,prior_matrix,nb_simul,use_seed=TRUE,seed_count=0, progressbarwidth=50)

    sd_simul=sapply(as.data.frame(rejection$summarystat), sd)
    
list(param=rejection$param, stats=rejection$summarystat, weights=array(1/nb_simul,nb_simul), stats_normalization=sd_simul, nsim=nb_simul, computime=as.numeric(difftime(Sys.time(), rejection$start, unit="secs")))
}


