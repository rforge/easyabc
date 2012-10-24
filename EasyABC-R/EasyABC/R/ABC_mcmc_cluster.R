## FUNCTION ABC_mcmc: ABC coupled to MCMC (Marjoram et al. 2003, Wegmann et al. 2009)
##############################################################################
ABC_mcmc_cluster <-function(method,model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,n_cluster=1,...){
    ## checking errors in the inputs
    if(missing(method)) stop("'method' is missing")
    if(missing(model)) stop("'model' is missing")
    if(missing(prior_matrix)) stop("'prior_matrix' is missing")
    if(missing(n_obs)) stop("'n_obs' is missing")
    if(missing(n_between_sampling)) stop("'n_between_sampling' is missing")
    if(missing(summary_stat_target)) stop("'summary_stat_target' is missing")
    if(!any(method == c("Marjoram_original", "Marjoram", "Wegmann"))){
        stop("Method must be Marjoram_original, Marjoram or wegmann")
    }
    if(!is.matrix(prior_matrix) && !is.data.frame(prior_matrix)) stop("'prior_matrix' has to be a matrix or data.frame.")
    if(is.data.frame(prior_matrix)) prior_matrix <- as.matrix(prior_matrix)
    if(dim(prior_matrix)[2]!=2) stop("'prior_matrix' must have two columns.")
    if(!is.vector(nb_simul)) stop("'nb_simul' has to be a number.")
    if(length(nb_simul)>1) stop("'nb_simul' has to be a number.")
    if (nb_simul<1) stop("'nb_simul' must be a number larger than 1.")
    nb_simul=floor(nb_simul)
    if(!is.vector(summary_stat_target)) stop("'summary_stat_target' has to be a vector.")
    if(!is.vector(n_cluster)) stop("'n_cluster' has to be a number.")
    if(length(n_cluster)>1) stop("'n_cluster' has to be a number.")
    if (n_cluster<1) stop ("'n_cluster' has to be a positive number.")
    n_cluster=floor(n_cluster)

	options(scipen=50)
	library(parallel)
## Note that we do not consider the original Marjoram's algortithm, which is not prone to parallel computing. (no calibration step)
	    return(switch(EXPR = method,
	       Marjoram = .ABC_MCMC2_cluster(model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,n_cluster,...),
	       Wegmann = .ABC_MCMC3_cluster(model,prior_matrix,n_obs,summary_stat_target,n_cluster,...)))

	options(scipen=0)
}
