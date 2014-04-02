library(EasyABC)
# library(lhs)
# library(mnormt)
# library(parallel)
# dyn.load("pkg/src/EasyABC.so")
# source("pkg/R/EasyABC-internal.R")
# source("pkg/R/ABC_sequential.R")
# source("pkg/R/trait_model.R")
# source("pkg/R/ABC_rejection.R")

test_explorations_equals <- function(ABC_sim1, ABC_sim2) {
  if (!all.equal(ABC_sim1$param,ABC_sim2$param)) {
    stop("Error in parameter sampling. They should be equals in the two explorations, with or without an additional constant.")
  }
  if (!all.equal(ABC_sim1$stats,ABC_sim2$stats)) {
    stop("Error in model evaluation or statistics computing. They should be equals in the two explorations, with or without an additional constant.")
  }
}

# compare the result with the original model and with an additional input variable fixed as constant
constants_test_rejection <- function(model1, prior1, model2, prior2, nb_simul) {
  set.seed(1)
  ABC_sim1<-ABC_rejection(model=model1, prior=prior1, nb_simul=nb_simul)
  set.seed(1)
  ABC_sim2<-ABC_rejection(model=model2, prior=prior2, nb_simul=nb_simul)
  test_explorations_equals(ABC_sim1, ABC_sim2)
  list(ABC_sim1=ABC_sim1,ABC_sim2=ABC_sim2)
}

# A model with a single argument
r = constants_test_rejection(
  model1=function(x){ 2 * x + 5 + rnorm(1,0,0.1) },
  prior1=list(c("unif",0,1)),
  model2=function(x){ x[2] * x[1] + 5 + rnorm(1,0,0.1) },
  prior2=list(c("unif",0,1),c("unif",2,2)),
  nb_simul=10)

# A model with two arguments
r = constants_test_rejection(
  model1=function(x){ c( x[1] + x[2] + rnorm(1,0,0.1) , x[1] * x[2] + rnorm(1,0,0.1) ) },
  prior1=list(c("unif",0,1),c("normal",1,2)),
  model2=function(x){ c( x[1] + x[2] * x[3] + rnorm(1,0,0.1) , x[1] * x[2] + rnorm(1,0,0.1) ) },
  prior2=list(c("unif",0,1),c("normal",1,2),c("unif",1,1)),
  nb_simul=10)
if (any(r$ABC_sim1$stats[,1]<0.5)||any(r$ABC_sim1$stats[,2]>4)) {
  stop("Wrong statistics with the test model")
}
  
constants_test_beaumont <- function(sum_stat_obs, tolerance_tab, model1, prior1, model2, prior2, nb_simul) {
  set.seed(1)
  ABC_sim1<-ABC_sequential(method="Beaumont", model=model1,
prior=prior1, nb_simul=nb_simul, summary_stat_target=sum_stat_obs,
tolerance_tab=tolerance_tab, use_seed=TRUE)
  set.seed(1)
  ABC_sim2<-ABC_sequential(method="Beaumont", model=model2,
prior=prior2, nb_simul=nb_simul, summary_stat_target=sum_stat_obs,
tolerance_tab=tolerance_tab, use_seed=TRUE)
  test_explorations_equals(ABC_sim1, ABC_sim2)
}

constants_test_beaumont(
  sum_stat_obs=c(100,2.5,20,30000),
  tolerance_tab=c(8,5),
  model1=trait_model,
  prior1=list(c("unif",3,5),c("unif",-2.3,1.6),c("unif",-25,125), c("unif",-0.7,3.2)),
  model2=function(input = c(1, 1, 1, 1, 1, 1)) {.C("trait_model", input = input, stat_to_return = array(0, 4))$stat_to_return} ,
  prior2=list(c("unif",500,500), c("unif",3,5),c("unif",-2.3,1.6), c("unif",1,1),c("unif",-25,125),c("unif",-0.7,3.2)),
  nb_simul=10)

