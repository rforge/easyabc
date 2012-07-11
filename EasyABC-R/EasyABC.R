######################################
## FUNCTIONS TO BE USED BY USERS (4)
######################################

## FUNCTION 1: brute-force ABC (Pritchard et al. 1999)
######################################################

ABC_rejection<-function(model,prior_matrix,nb_simul,use_seed=TRUE,seed_count=0){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	for (i in 1:nb_simul){
		l=dim(prior_matrix)[1]
		param=array(0,l)
		for (j in 1:l){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		if (use_seed) {
			param=c((seed_count+i),param)
		}
		simul_summarystat=model(param)
		tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
		if (use_seed) {
			tab_param=rbind(tab_param,param[2:(l+1)])
		}
		else{
			tab_param=rbind(tab_param,param)
		}
	}
	cbind(tab_param,tab_simul_summarystat)
}

## FUNCTION 2: Sequential ABC methods (Beaumont et al. 2009, Drovandi & Pettitt 2011, Del Moral et al. 2011, Lenormand et al. 2012)
###################################################################################################################################
ABC_sequential<-function(){
	## fonction à faire qui regroupe les algos sequentiels (PMC, drovandi, delmoral, maxime)
}

## FUNCTION 3: ABC coupled to MCMC (Marjoram et al. 2003, Wegmann et al. 2009)
##############################################################################
ABC_mcmc<-function(){
	## fonction à faire qui regroupe les algos de mcmc (marjoram, marjoram modifié, wegmann)
}

## FUNCTION 4: ABC coupled to emulators (XXX)
#############################################
ABC_emulator<-function(){
	## fonction à faire qui regroupe les algos de type emulator
}

## sample usage
###############
prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
prior_matrix
# linux
# ABC_rejection(.binary_model("./parthy"),prior_matrix,10,TRUE)
# windows
# ABC_rejection(.binary_model("./parthy_test.exe"),prior_matrix,10,TRUE)


#######################
## INTERNAL FUNCTIONS
#######################

## model wrapper
################
.binary_model<-function(command) {
	invoke<-function(param) {
		write.table(param,file="input",row.names=F,col.names=F,quote=F)
		system(command)
		read.table("output",h=F)
	}
}

## function to compute a distance between a matrix of simulated statistics (row: different simulations, columns: different summary statistics) and the array of data summary statistics
#######################################################################################################################################################################################
.compute_dist<-function(summary_stat_target,simul,sd_simul){
	l=length(summary_stat_target)
	nsimul=dim(simul)[1]
	vartab=array(1,l)
	dist=array(0,nsimul)
	for (i in 1:l){
		vartab[i]=min(1,1/(sd_simul[i]*sd_simul[i])) ## differences between simul and data are normalized in each dimension by the empirical variances in each dimension
		dist=dist+vartab[i]*(simul[,i]-summary_stat_target[i])*(simul[,i]-summary_stat_target[i]) ## an euclidean distance is used
	}
dist
}

## same as .compute_dist when there is only one simulation
##########################################################
.compute_dist_single<-function(summary_stat_target,simul,sd_simul){
	l=length(summary_stat_target)
	dist=0
	vartab=array(1,l)
	for (i in 1:l){
		vartab[i]=min(1,1/(sd_simul[i]*sd_simul[i]))
		dist=dist+vartab[i]*(simul[i]-summary_stat_target[i])*(simul[i]-summary_stat_target[i])
	}
dist
}

## function to select the simulations that are at a distance smaller than tol from the data
###########################################################################################
.selec_simul<-function(summary_stat_target,param,simul,sd_simul,tol){
	dist=.compute_dist(summary_stat_target,simul,sd_simul)
	res=cbind(param[dist<tol,],simul[dist<tol,])
	if (dim(res)[1]==0){
		res=NULL
	}
res
}

## function to randomly pick a particle from a weighted array (of sum=1)
########################################################################
.particle_pick<-function(param,tab_weight){
	u=runif(1)
	weight_cum=cumsum(tab_weight)
	pos=1:length(tab_weight)
	p=min(pos[weight_cum>u])
param[p,]
}

# with package mnormt
library(mnormt)

## function to move a particle
##############################
.move_particle<-function(param_picked,varcov_matrix){
rmnorm(n = 1, mean = param_picked, varcov_matrix)
}

## function to perform ABC simulations from a non-uniform prior (derived from a set of particles)
#################################################################################################
.ABC_launcher_not_uniform<-function(model,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,use_seed,seed_count){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	for (i in 1:nb_simul){
		l=dim(param_previous_step)[2]
		# pick a particle
		param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
		# move it
		param_moved=.move_particle(param_picked,2*var(param_previous_step[,tab_unfixed_param])) # only variable parameters are moved
		param=param_previous_step[1,]
		param[tab_unfixed_param]=param_moved
		if (use_seed) {
			param=c((seed_count+i),param)
		}
		simul_summarystat=model(param)
		tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
		if (use_seed) {
			tab_param=rbind(tab_param,param[2:(l+1)])
		}
		else{
			tab_param=rbind(tab_param,param)
		}
	}
	cbind(tab_param,tab_simul_summarystat)
}

## function to compute particle weights
####################################### 
.compute_weight<-function(param_simulated,param_previous_step,tab_weight){
	vmat=2*var(param_previous_step)
	n_particle=dim(param_previous_step)[1]
	n_new_particle=dim(param_simulated)[1]
	tab_weight_new=array(0,n_new_particle)
	for (i in 1:n_particle){
		for (j in 1:n_new_particle){
			tab_weight_new[j]=tab_weight_new[j]+tab_weight[i]*dmnorm(param_simulated[j,],param_previous_step[i,],vmat)
		}
	}
	tab_weight_new=1/tab_weight_new
tab_weight_new/sum(tab_weight_new)
}

## PMC ABC algorithm: Beaumont et al. Biometrika 2009
#####################################################
.ABC_PMC<-function(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed=TRUE){
	T=length(tolerance_tab)
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)
	seed_count=0
	tab_unfixed_param=array(TRUE,nparam)
	for (i in 1:nparam){
		tab_unfixed_param[i]=(prior_matrix[i,1]!=prior_matrix[i,2])
	}

## step 1
	nb_simul_step=nb_simul
	simul_below_tol=NULL
	while (nb_simul_step>0){
		if (nb_simul_step>1){
			# classic ABC step
			tab_ini=ABC_rejection(model,prior_matrix,nb_simul_step,use_seed,seed_count)
			sd_simul=sd(tab_ini[,(nparam+1):(nparam+nstat)]) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
			seed_count=seed_count+nb_simul_step
			# selection of simulations below the first tolerance level
			simul_below_tol=rbind(simul_below_tol,.selec_simul(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,tolerance_tab[1]))
			if (length(simul_below_tol)>0){
				nb_simul_step=nb_simul-dim(simul_below_tol)[1]
			}
		}
		else{
			tab_ini=ABC_rejection(model,prior_matrix,nb_simul_step,use_seed,seed_count)
			seed_count=seed_count+nb_simul_step
			if (.compute_dist_single(summary_stat_target,tab_ini[(nparam+1):(nparam+nstat)],sd_simul)<tolerance_tab[1]){
				simul_below_tol=rbind(simul_below_tol,tab_ini)
				nb_simul_step=0
			}
		}
	} # until we get nb_simul simulations below the first tolerance threshold
	# initially, weights are equal
	tab_weight=array(1/nb_simul,nb_simul)
	write.table(cbind(tab_weight,simul_below_tol),file="output_step1",row.names=F,col.names=F,quote=F)

## steps 2 to T
	for (it in 2:T){
		nb_simul_step=nb_simul
		simul_below_tol2=NULL
		while (nb_simul_step>0){
			if (nb_simul_step>1){
				# Sampling of parameters around the previous particles
				tab_ini=.ABC_launcher_not_uniform(model,simul_below_tol[,1:nparam],tab_unfixed_param,tab_weight,nb_simul_step,use_seed,seed_count)
				seed_count=seed_count+nb_simul_step
				simul_below_tol2=rbind(simul_below_tol2,.selec_simul(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,tolerance_tab[it]))
				if (length(simul_below_tol2)>0){
					nb_simul_step=nb_simul-dim(simul_below_tol2)[1]
				}
			}
			else{
				tab_ini=ABC_rejection(model,prior_matrix,nb_simul_step,use_seed,seed_count)
				seed_count=seed_count+nb_simul_step
				if (.compute_dist_single(summary_stat_target,tab_ini[(nparam+1):(nparam+nstat)],sd_simul)<tolerance_tab[it]){
					simul_below_tol=rbind(simul_below_tol,tab_ini)
					nb_simul_step=0
				}
			}
		} # until we get nb_simul simulations below the it-th tolerance threshold
		# update of particle weights
		tab_weight2=.compute_weight(simul_below_tol2[,1:nparam],simul_below_tol[,1:nparam],tab_weight)
		# update of the set of particles and of the associated weights for the next ABC sequence
		tab_weight=tab_weight2
		simul_below_tol=simul_below_tol2
		write.table(cbind(tab_weight,simul_below_tol),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)

	}
cbind(tab_weight,simul_below_tol)
}

## test
# linux
# .ABC_PMC(binary_model("./parthy"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5),use_seed=TRUE)
# windows
# .ABC_PMC(binary_model("./parthy_test.exe"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5),use_seed=TRUE)


## Algo de Marjoram

## Algo de Marjoram avec détermination de epsilon et proposal_range de Wegmann et al. 2009
.ABC_MCMC2<-function(model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,proposal_method="unif",burn_in_length=0,n_calibration=10000,tolerance_quantile=0.01,proposal_phi=1,use_seed=FALSE){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	# initial draw of a particle
	for (i in 1:(n_calibration)){
		l=dim(prior_matrix)[1]
		param=array(0,l)
		for (j in 1:l){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		if (use_seed) {
			param=c(i,param)
		}
		simul_summarystat=model(param)
		tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
		tab_param=rbind(tab_param,param)
	}
	simuldist=.compute_dist(summary_stat_target,tab_simul_summarystat,sd_simul)
	ord_sim=order(simuldist,decreasing=F)
	nmax=ceiling(tolerance_quantile*n_calibration)
	dist_max=simuldist[(ord_sim[nmax])]
	proposal_range=array(0,l)
	tab_param=tab_param[(ord_sim[1:nmax]),]
	n_ini=sample(nmax,1)
	param_ini=tab_param[n_ini,]
	tab_simul_ini=tab_simul_summarystat[(ord_sim[n_ini]),]
	if (use_seed) {
			tab_param=tab_param[,2:(l+1)]
	}
	for (i in 1:l){
		proposal_range[l]=sd(tab_param[,l])*proposal_phi
	}

	# chain run
	tab_param=param_ini
	tab_simul_summarystat=tab_simul_ini
	seed_compt=n_calibration+1
	for (is in 2:n_obs){
		for (i in 1:n_between_sampling){
			param=.move_particle(param_ini,proposal_method,proposal_range,use_seed,seed_compt)
			simul_summarystat=model(param)
			dist_simul=.compute_dist_single(summary_stat_target,simul_summary_stat,proposal_range/phi)
			if (dist_simul<dist_max){
				param_ini=param
				tab_simul_ini=simul_summarystat
			}
			seed_compt=seed_compt+1
		}
		tab_simul_summarystat=rbind(tab_simul_summarystat,tab_simul_ini)
		tab_param=rbind(tab_param,param_ini)
	}	
cbind(tab_param,tab_simul_summarystat)	
}

