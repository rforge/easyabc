## FUNCTION ABC_mcmc: ABC coupled to MCMC (Marjoram et al. 2003, Wegmann et al.
## 2009)
ABC_mcmc <- function(method, model, prior, summary_stat_target, prior_test = NULL, 
    n_rec = 100, n_between_sampling = 10, n_cluster = 1, use_seed = FALSE, verbose = FALSE, 
    ...) {
    ## checking errors in the inputs
    if (missing(method)) 
        stop("'method' is missing")
    if (missing(model)) 
        stop("'model' is missing")
    if (missing(prior)) 
        stop("'prior' is missing")
    if (!is.null(prior_test)) 
        .check_prior_test(length(prior), prior_test)
    data = .wrap_constants_in_model(prior, model, use_seed)
    prior = data$new_prior
    model = data$new_model
    prior = .process_prior(prior)
    if (!any(method == c("Marjoram_original", "Marjoram", "Wegmann"))) {
        stop("Method must be Marjoram_original, Marjoram or wegmann")
    }
    .test_param("summary_stat_target", missing = TRUE, vector = TRUE)
    .test_param("n_cluster", positive = TRUE)
    n_cluster = floor(n_cluster)
    if (!is.logical(use_seed)) 
        stop("'use_seed' has to be boolean")
    if (!is.logical(verbose)) 
        stop("'verbose' has to be boolean")
    mcmc = NULL
    if (n_cluster == 1) {
        mcmc = .ABC_mcmc_internal(method, model, prior, prior_test, n_rec, n_between_sampling, 
            summary_stat_target, use_seed, verbose, ...)
    } else {
        mcmc = .ABC_mcmc_cluster(method, model, prior, prior_test, n_rec, n_between_sampling, 
            summary_stat_target, n_cluster, use_seed, verbose, ...)
    }
    mcmc
} 
