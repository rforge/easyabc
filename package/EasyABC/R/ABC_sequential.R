ABC_sequential <-
function(functionName,model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed=TRUE,seed_count=0,inside_prior=TRUE){
	## fonction regroup sequentials algorithms (PMC, drovandi, delmoral, maxime)
  switch(EXPR = functionName, 
         PMC = .ABC_PMC(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed,seed_count,inside_prior),
         PMC2 = .ABC_PMC2(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed,seed_count,inside_prior),
         Drovandi = .ABC_Drovandi(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,alpha,c,first_tolerance_level_auto,use_seed,seed_count),
         Delmoral = .ABC_Delmoral(model,prior_matrix,nb_simul,alpha,M,nb_threshold,tolerance_target,summary_stat_target,use_seed,seed_count),
         Lenormand = .ABC_Lenormand(model,prior_matrix,nb_simul,summary_stat_target,alpha,p_acc_min,use_seed,seed_count,inside_prior))
}

