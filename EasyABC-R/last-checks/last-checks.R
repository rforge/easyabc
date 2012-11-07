library(EasyABC)
set.seed(1)
priormatrix=cbind(c(500,3,-2.3,1,-25,-0.7),c(500,5,1.6,1,125,3.2))
priormatrix
sum_stat_obs=c(100,2.5,20,30000)
n=1000

### Rejection
##############

## Rejection - simple core
set.seed(1)
ABC_rej<-ABC_rejection(model=trait_model, prior_matrix=priormatrix, nb_simul=n)
ABC_rej$computime
write.table(cbind(ABC_rej$weights,ABC_rej$param,ABC_rej$stats),file="ABC_rejection_simple_core",col.names=F,row.names=F,quote=F)

## Rejection - multiple cores
set.seed(1)
ABC_rejb<-ABC_rejection(model=trait_model, prior_matrix=priormatrix, nb_simul=n,n_cluster=2)
ABC_rejb$computime
write.table(cbind(ABC_rejb$weights,ABC_rejb$param,ABC_rejb$stats),file="ABC_rejection_multiple_cores",col.names=F,row.names=F,quote=F)



##################
### sequential ###
##################

### Beaumont
#############

## Beaumont - simple core
set.seed(1)
tolerance=c(2.5,2,1.75)
ABC_Beaumont<-ABC_sequential(method="Beaumont", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs,tolerance_tab=tolerance)
ABC_Beaumont$computime
write.table(cbind(ABC_Beaumont$weights,ABC_Beaumont$param,ABC_Beaumont$stats),file="ABC_Beaumont_simple_core",col.names=F,row.names=F,quote=F)

## Beaumont - multiple cores
set.seed(1)
tolerance=c(2.5,2,1.75)
ABC_Beaumontb<-ABC_sequential(method="Beaumont", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs,tolerance_tab=tolerance,n_cluster=2)
ABC_Beaumontb$computime
write.table(cbind(ABC_Beaumontb$weights,ABC_Beaumontb$param,ABC_Beaumontb$stats),file="ABC_Beaumont_multiple_cores",col.names=F,row.names=F,quote=F)



### Drovandi
#############

## Drovandi - simple core
set.seed(1)
tolerance=1.75
ABC_Drovandi<-ABC_sequential(method="Drovandi", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs,tolerance_tab=tolerance, c=0.05)
ABC_Drovandi$computime
write.table(cbind(ABC_Drovandi$weights,ABC_Drovandi$param,ABC_Drovandi$stats),file="ABC_Drovandi_simple_core",col.names=F,row.names=F,quote=F)

## Drovandi - multiple cores
set.seed(1)
tolerance=1.75
ABC_Drovandib<-ABC_sequential(method="Drovandi", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs,tolerance_tab=tolerance, c=0.05, n_cluster=2)
ABC_Drovandib$computime
write.table(cbind(ABC_Drovandib$weights,ABC_Drovandib$param,ABC_Drovandib$stats),file="ABC_Drovandi_multiple_cores",col.names=F,row.names=F,quote=F)



### Del Moral
##############

## Del Moral - simple core
set.seed(1)
tolerance=2.25
ABC_Delmoral<-ABC_sequential(method="Delmoral", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs,tolerance_target=tolerance)
ABC_Delmoral$computime
write.table(cbind(ABC_Delmoral$weights,ABC_Delmoral$param,ABC_Delmoral$stats),file="ABC_Delmoral_simple_core",col.names=F,row.names=F,quote=F)

## Delmoral - multiple cores
set.seed(1)
tolerance=2.25
ABC_Delmoralb<-ABC_sequential(method="Delmoral", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs,tolerance_target=tolerance, n_cluster=2)
ABC_Delmoralb$computime
write.table(cbind(ABC_Delmoralb$weights,ABC_Delmoralb$param,ABC_Delmoralb$stats),file="ABC_Delmoral_multiple_cores",col.names=F,row.names=F,quote=F)

## Del Moral - M=15 - simple core
set.seed(1)
tolerance=3.5
ABC_Delmoral15<-ABC_sequential(method="Delmoral", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs,tolerance_target=tolerance,M=15)
ABC_Delmoral15$computime
write.table(cbind(ABC_Delmoral15$weights,ABC_Delmoral15$param,ABC_Delmoral15$stats),file="ABC_Delmoral_M15_simple_core",col.names=F,row.names=F,quote=F)

## Delmoral - M=15 - multiple cores
set.seed(1)
tolerance=3.5
ABC_Delmoral15b<-ABC_sequential(method="Delmoral", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs,tolerance_target=tolerance, M=15, n_cluster=2)
ABC_Delmoral15b$computime
write.table(cbind(ABC_Delmoral15b$weights,ABC_Delmoral15b$param,ABC_Delmoral15b$stats),file="ABC_Delmoral_M15_multiple_cores",col.names=F,row.names=F,quote=F)



### Lenormand
##############

## Lenormand - simple core
set.seed(1)
paccmin=0.1
ABC_Lenormand<-ABC_sequential(method="Lenormand", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs, p_acc_min=paccmin)
ABC_Lenormand$computime
write.table(cbind(ABC_Lenormand$weights,ABC_Lenormand$param,ABC_Lenormand$stats),file="ABC_Lenormand_simple_core",col.names=F,row.names=F,quote=F)

## Lenormand - multiple cores
set.seed(1)
paccmin=0.1
ABC_Lenormandb<-ABC_sequential(method="Lenormand", model=trait_model,prior_matrix=priormatrix, nb_simul=n, summary_stat_target=sum_stat_obs, p_acc_min=paccmin, n_cluster=2)
ABC_Lenormandb$computime
write.table(cbind(ABC_Lenormandb$weights,ABC_Lenormandb$param,ABC_Lenormandb$stats),file="ABC_Lenormand_multiple_cores",col.names=F,row.names=F,quote=F)



############
### MCMC ###
############

### Marjoram original
#############

## Marjoram original- simple core
set.seed(1)
nbetweensampling=1
distmax=2.5
tabnormalization=c(50,1,20,10000)
proposalrange=c(0,1,0.5,0,50,1)
ABC_Marjoram_original<-ABC_mcmc(method="Marjoram_original",model=trait_model,prior_matrix=priormatrix,n_obs=n,n_between_sampling=nbetweensampling,summary_stat_target=sum_stat_obs,dist_max=distmax,tab_normalization=tabnormalization,proposal_range=proposalrange)
ABC_Marjoram_original$computime
write.table(cbind(ABC_Marjoram_original$param,ABC_Marjoram_original$stats),file="ABC_Marjoram_original_simple_core",col.names=F,row.names=F,quote=F)

## Marjoram original- multiple cores
set.seed(1)
nbetweensampling=1
distmax=2.5
tabnormalization=c(50,1,20,10000)
proposalrange=c(0,1,0.5,0,50,1)
ABC_Marjoram_originalb<-ABC_mcmc(method="Marjoram_original",model=trait_model,prior_matrix=priormatrix,n_obs=n,n_between_sampling=nbetweensampling,summary_stat_target=sum_stat_obs,dist_max=distmax,tab_normalization=tabnormalization,proposal_range=proposalrange,n_cluster=2)
ABC_Marjoram_originalb$computime
write.table(cbind(ABC_Marjoram_originalb$param,ABC_Marjoram_originalb$stats),file="ABC_Marjoram_original_multiple_cores",col.names=F,row.names=F,quote=F)


### Marjoram
#############

## Marjoram - simple core
set.seed(1)
nbetweensampling=1
ncalib=1000
tolquantile = 0.05
proposalphi=1
ABC_Marjoram<-ABC_mcmc(method="Marjoram",model=trait_model,prior_matrix=priormatrix,n_obs=n,n_between_sampling=nbetweensampling,summary_stat_target=sum_stat_obs,n_calibration=ncalib,tolerance_quantile=tolquantile,proposal_phi=proposalphi)
ABC_Marjoram$computime
write.table(cbind(ABC_Marjoram$param,ABC_Marjoram$stats),file="ABC_Marjoram_simple_core",col.names=F,row.names=F,quote=F)

## Marjoram - multiple cores
set.seed(1)
nbetweensampling=1
ncalib=1000
tolquantile = 0.05
proposalphi=1
ABC_Marjoramb<-ABC_mcmc(method="Marjoram",model=trait_model,prior_matrix=priormatrix,n_obs=n,n_between_sampling=nbetweensampling,summary_stat_target=sum_stat_obs,n_calibration=ncalib,tolerance_quantile=tolquantile,proposal_phi=proposalphi,n_cluster=2)
ABC_Marjoramb$computime
write.table(cbind(ABC_Marjoramb$param,ABC_Marjoramb$stats),file="ABC_Marjoram_multiple_cores",col.names=F,row.names=F,quote=F)


### Wegmann
############

## Wegmann - simple core
set.seed(1)
nbetweensampling=1
ncalib=1000
tolquantile = 0.05
proposalphi=1
ABC_Wegmann<-ABC_mcmc(method="Wegmann",model=trait_model,prior_matrix=priormatrix,n_obs=n,n_between_sampling=nbetweensampling,summary_stat_target=sum_stat_obs,n_calibration=ncalib,tolerance_quantile=tolquantile,proposal_phi=proposalphi,numcomp=0)
ABC_Wegmann$computime
write.table(cbind(ABC_Wegmann$param,ABC_Wegmann$stats),file="ABC_Wegmann_simple_core",col.names=F,row.names=F,quote=F)

## Wegmann - multiple cores
set.seed(1)
nbetweensampling=1
ncalib=1000
tolquantile = 0.05
proposalphi=1
ABC_Wegmannb<-ABC_mcmc(method="Wegmann",model=trait_model,prior_matrix=priormatrix,n_obs=n,n_between_sampling=nbetweensampling,summary_stat_target=sum_stat_obs,n_calibration=ncalib,tolerance_quantile=tolquantile,proposal_phi=proposalphi,numcomp=0,n_cluster=2)
ABC_Wegmannb$computime
write.table(cbind(ABC_Wegmannb$param,ABC_Wegmannb$stats),file="ABC_Wegmann_multiple_cores",col.names=F,row.names=F,quote=F)

