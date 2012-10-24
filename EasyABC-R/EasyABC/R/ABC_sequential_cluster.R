ABC_sequential_cluster <-
function(method,model,prior_matrix,nb_simul,summary_stat_target,n_cluster=1,...){
    ## checking errors in the inputs
    if(missing(method)) stop("'method' is missing")
    if(missing(model)) stop("'model' is missing")
    if(missing(prior_matrix)) stop("'prior_matrix' is missing")
    if(missing(nb_simul)) stop("'nb_simul' is missing")
    if(missing(summary_stat_target)) stop("'summary_stat_target' is missing")
    if(!any(method == c("Beaumont", "Drovandi", "Delmoral","Lenormand"))){
        stop("Method must be Beaumont, Drovandi, Delmoral or Lenormand")
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

	       ## general function regrouping the different sequential algorithms     
	  ## [Beaumont et al., 2009] Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009). Adaptive approximate Bayesian computation. Biometrika,96(4):983â??990.
	  ## [Drovandi & Pettitt 2011] Drovandi, C. C. and Pettitt, A. N. (2011). Estimation of parameters for macroparasite population evolution using approximate Bayesian computation. Biometrics, 67(1):225â??233.
	  ## [Del Moral et al. 2012] Del Moral, P., Doucet, A., and Jasra, A. (2012). An adaptive sequential Monte Carlo method for approximate Bayesian computation, Statistics and Computing., 22(5):1009-1020.
	  ## [Lenormand et al. 2012] Lenormand, M., Jabot, F., Deffuant G. (2012). Adaptive approximate Bayesian computation for complex models, submitted to Comput. Stat. )
	  switch(EXPR = method,
         	 Beaumont = .ABC_PMC_cluster(model,prior_matrix,nb_simul,summary_stat_target,n_cluster,,...),
	         Drovandi = .ABC_Drovandi_cluster(model,prior_matrix,nb_simul,summary_stat_target,n_cluster,...),
	         Delmoral = .ABC_Delmoral_cluster(model,prior_matrix,nb_simul,summary_stat_target,n_cluster,...),
	         Lenormand = .ABC_Lenormand_cluster(model,prior_matrix,nb_simul,summary_stat_target,n_cluster,...))

	options(scipen=0)

}
