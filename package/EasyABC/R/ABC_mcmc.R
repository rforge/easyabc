ABC_mcmc <-
	## fonction ? faire qui regroupe les algos de mcmc (marjoram, marjoram modifi?, wegmann)
  function(functionName,model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,dist_max,tab_normalization,proposal_range,use_seed=TRUE,seed_count=0,...){
    switch(EXPR = functionName, 
       MCMC = .ABC_MCMC(model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,dist_max,tab_normalization,proposal_range,use_seed,seed_count),
       MCMC2 = .ABC_MCMC2(model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,n_calibration,tolerance_quantile,proposal_phi,use_seed,seed_count),
       MCMC3 = .ABC_MCMC3(model,prior_matrix,n_obs,summary_stat_targ,n_calibration,tolerance_quantile,proposal_phi,n_between_sampling,numcomp,use_seed,seed_count))
  }


