ABC_sequential <-
function(functionName,model,prior_matrix,nb_simul,summary_stat_target,...){
	## fonction regroup sequentials algorithms     
  ## [Beaumont et al., 2009] Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009). Adaptive approximate Bayesian computation. Biometrika,96(4):983–990.
  ## [Drovandi & Pettitt 2011] Drovandi, C. C. and Pettitt, A. N. (2011). Estimation of parameters for macroparasite population evolution using approximate Bayesian computation. Biometrics, 67(1):225–233.
  ## [Del Moral et al. 2012] Del Moral, P., Doucet, A., and Jasra, A. (2012). An adaptive sequential Monte Carlo method for approximate Bayesian computation, Statistics and Computing., 22(5):1009-1020.
  ## [Lenormand et al. 2012] Lenormand, M., Jabot, F., Deffuant G. (2012). Adaptive approximate Bayesian computation for complex models, submitted to Comput. Stat. )
  switch(EXPR = functionName, 
         PMC = .ABC_PMC(model,prior_matrix,nb_simul,summary_stat_target,...),
         PMC2 = .ABC_PMC2(model,prior_matrix,nb_simul,summary_stat_target,...),
         Drovandi = .ABC_Drovandi(model,prior_matrix,nb_simul,summary_stat_target,...),
         Delmoral = .ABC_Delmoral(model,prior_matrix,nb_simul,summary_stat_target,...),
         Lenormand = .ABC_Lenormand(model,prior_matrix,nb_simul,summary_stat_target,...))
}

