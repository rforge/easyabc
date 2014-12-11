## FUNCTION ABC_sequential: Sequential ABC methods (Beaumont et al. 2009, Drovandi
## & Pettitt 2011, Del Moral et al. 2011, Lenormand et al. 2012)
ABC_sequential <- function(method, model, prior, nb_simul, summary_stat_target, prior_test = NULL, 
    n_cluster = 1, use_seed = FALSE, verbose = FALSE, ...) {
    ## checking errors in the inputs
    if (missing(method)) 
        stop("'method' is missing")
    if (missing(model)) 
        stop("'model' is missing")
    if (missing(prior)) 
        stop("'prior' is missing")
    data = .wrap_constants_in_model(prior, model, use_seed)
    prior = data$new_prior
    model = data$new_model
    prior = .process_prior(prior)
    if (!is.null(prior_test)) 
        .check_prior_test(length(prior), prior_test)
    .test_param("summary_stat_target", missing = TRUE, vector = TRUE)
    if (!any(method == c("Beaumont", "Drovandi", "Delmoral", "Lenormand", "Emulation"))) {
        stop("Method must be Beaumont, Drovandi, Delmoral, Lenormand or Emulation")
    }
    .test_param("nb_simul", missing = TRUE)
    if (nb_simul < 1) 
        stop("'nb_simul' must be a number larger than 1.")
    nb_simul = floor(nb_simul)
    .test_param("n_cluster", positive = TRUE)
    n_cluster = floor(n_cluster)
    if (!is.logical(use_seed)) 
        stop("'use_seed' has to be boolean")
    if (!is.logical(verbose)) 
        stop("'verbose' has to be boolean")
    sequential = NULL
    if (n_cluster == 1) {
        sequential = .ABC_sequential(method, model, prior, prior_test, nb_simul, 
            summary_stat_target, use_seed, verbose, ...)
    } else {
        sequential = .ABC_sequential_cluster(method, model, prior, prior_test, nb_simul, 
            summary_stat_target, n_cluster, use_seed, verbose, ...)
    }
    sequential
} 
