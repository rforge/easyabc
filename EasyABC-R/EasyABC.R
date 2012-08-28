######################################
## FUNCTIONS TO BE USED BY USERS (4)
######################################

## FUNCTION 1: brute-force ABC (Pritchard et al. 1999)
######################################################

ABC_rejection<-function(model,prior_matrix,nb_simul,use_seed=TRUE,seed_count=0){
	options(scipen=50)
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
	options(scipen=0)
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


prior_matrix=c(500,3,-2.3,1,-0.25,-0.7,500,5,1.6,1,1.25,3.2)
dim(prior_matrix)<-c(6,2)
prior_matrix
# linux
# ABC_rejection(.binary_model("./trait_model"),prior_matrix,10,TRUE)
# windows
# ABC_rejection(.binary_model("./trait_model.exe"),prior_matrix,10,TRUE)

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

## function to check whether moved parameters are still in the prior distribution
#################################################################################
.is_included<-function(res,prior_matrix){
	test=TRUE
	for (i in 1:length(res)){
		if ((res[i]<prior_matrix[i,1])||(res[i]>prior_matrix[i,2])){
			test=FALSE
		}
	}
test	
}

## function to move a particle
##############################
.move_particle<-function(param_picked,varcov_matrix,prior_matrix){
	rmnorm(n = 1, mean = param_picked, varcov_matrix)
}

## function to move a particle
##############################
.move_particleb<-function(param_picked,varcov_matrix,prior_matrix){
	test=FALSE
	while (!test){
		res=rmnorm(n = 1, mean = param_picked, varcov_matrix)
		test=.is_included(res,prior_matrix)
	}
res
}


## function to perform ABC simulations from a non-uniform prior (derived from a set of particles)
#################################################################################################
.ABC_launcher_not_uniform<-function(model,prior_matrix,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,use_seed,seed_count,inside_prior){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	for (i in 1:nb_simul){
		l=dim(param_previous_step)[2]
		if (!inside_prior){
			# pick a particle
			param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
			# move it
			param_moved=.move_particle(param_picked,2*cov.wt(param_previous_step[,tab_unfixed_param],as.vector(tab_weight))$cov,prior_matrix[tab_unfixed_param,]) # only variable parameters are moved, computation of a WEIGHTED variance
		}
		else{
			test=FALSE
			while (!test){
				# pick a particle
				param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
				# move it
				param_moved=.move_particle(param_picked,2*cov.wt(param_previous_step[,tab_unfixed_param],as.vector(tab_weight))$cov,prior_matrix[tab_unfixed_param,]) # only variable parameters are moved, computation of a WEIGHTED variance
				test=.is_included(param_moved,prior_matrix)
			}

		}
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

## PMC ABC algorithm with multivariate normal jumps
###################################################
.ABC_PMC2<-function(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed=TRUE,seed_count=0,inside_prior=TRUE){
	T=length(tolerance_tab)
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)
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
	print("step 1 completed")
## steps 2 to T
	for (it in 2:T){
		nb_simul_step=nb_simul
		simul_below_tol2=NULL
		while (nb_simul_step>0){
			if (nb_simul_step>1){
				# Sampling of parameters around the previous particles
				tab_ini=.ABC_launcher_not_uniform(model,prior_matrix,simul_below_tol[,1:nparam],tab_unfixed_param,tab_weight,nb_simul_step,use_seed,seed_count,inside_prior)
				seed_count=seed_count+nb_simul_step
				simul_below_tol2=rbind(simul_below_tol2,.selec_simul(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,tolerance_tab[it]))
				if (length(simul_below_tol2)>0){
					nb_simul_step=nb_simul-dim(simul_below_tol2)[1]
				}
			}
			else{
				tab_ini=.ABC_launcher_not_uniform(model,prior_matrix,simul_below_tol[,1:nparam],tab_unfixed_param,tab_weight,nb_simul_step,use_seed,seed_count,inside_prior)
				seed_count=seed_count+nb_simul_step
				if (.compute_dist_single(summary_stat_target,tab_ini[(nparam+1):(nparam+nstat)],sd_simul)<tolerance_tab[it]){
					simul_below_tol2=rbind(simul_below_tol2,tab_ini)
					nb_simul_step=0
				}
			}
		} # until we get nb_simul simulations below the it-th tolerance threshold
		# update of particle weights
		tab_weight2=.compute_weight(simul_below_tol2[,1:nparam][,tab_unfixed_param],simul_below_tol[,1:nparam][,tab_unfixed_param],tab_weight)
		# update of the set of particles and of the associated weights for the next ABC sequence
		tab_weight=tab_weight2
		simul_below_tol=matrix(0,nb_simul,(nparam+nstat))
		for (i1 in 1:nb_simul){
			for (i2 in 1:(nparam+nstat)){
				simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
			}
		}
		write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
		print(paste("step ",it," completed",sep=""))
	}
cbind(tab_weight,simul_below_tol)
}

## test
# linux
# .ABC_PMC2(.binary_model("./parthy"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5),use_seed=TRUE,inside_prior=TRUE)
# windows
# .ABC_PMC2(.binary_model("./parthy_test.exe"),prior_matrix,20,c(1,0.9,0.8),c(50,2.5),use_seed=TRUE,inside_prior=TRUE)

## function to move a particle with a unidimensional normal jump
################################################################
.move_particle_uni<-function(param_picked,sd_array,prior_matrix){
	res=param_picked
	for (i in 1:length(param_picked)){
		res[i]=rnorm(n = 1, mean = param_picked[i], sd_array[i])
	}
res
}

## function to move a particle with a unidimensional normal jump
################################################################
.move_particleb_uni<-function(param_picked,sd_array,prior_matrix){
	test=FALSE
	res=param_picked
	while (!test){
		for (i in 1:length(param_picked)){
			res[i]=rnorm(n = 1, mean = param_picked[i], sd_array[i])
		}
		test=.is_included(res,prior_matrix)
	}
res
}


## function to compute particle weights with unidimensional jumps
#################################################################
.compute_weight_uni<-function(param_simulated,param_previous_step,tab_weight){
	l=dim(param_previous_step)[2]
	sd_array=array(1,l)
	for (j in 1:l){
		sd_array[j]=sqrt(2*var(param_previous_step[,j]))
	}
	sd_array=as.numeric(sd_array)
	n_particle=dim(param_previous_step)[1]
	n_new_particle=dim(param_simulated)[1]
	tab_weight_new=array(0,n_new_particle)
	for (i in 1:n_particle){
		for (j in 1:n_new_particle){
			tab_temp=tab_weight[i]
			for (k in 1:l){
				tab_temp=tab_temp*dnorm(as.numeric(param_simulated[j,k]),as.numeric(param_previous_step[i,k]),sd_array[k])
			}
			tab_weight_new[j]=tab_weight_new[j]+tab_temp
		}
	}
	tab_weight_new=1/tab_weight_new
tab_weight_new/sum(tab_weight_new)
}

## function to perform ABC simulations from a non-uniform prior and with unidimensional jumps
#############################################################################################
.ABC_launcher_not_uniform_uni<-function(model,prior_matrix,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,use_seed,seed_count,inside_prior){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	for (i in 1:nb_simul){
		l=dim(param_previous_step)[2]
		l_array=dim(param_previous_step[,tab_unfixed_param])[2]
		sd_array=array(1,l_array)
		covmat=2*cov.wt(param_previous_step[,tab_unfixed_param],as.vector(tab_weight))$cov # computation of a WEIGHTED variance
		for (j in 1:l_array){
			sd_array[j]=sqrt(covmat[j,j])

		}
		if (!inside_prior){
			# pick a particle
			param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
			# move it
			param_moved=.move_particle_uni(param_picked,sd_array,prior_matrix[tab_unfixed_param,]) # only variable parameters are moved
		}
		else{
			test=FALSE
			while (!test){
				# pick a particle
				param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
				# move it
				param_moved=.move_particle_uni(param_picked,sd_array,prior_matrix[tab_unfixed_param,]) # only variable parameters are moved
				test=.is_included(param_moved,prior_matrix)
			}

		}
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

## PMC ABC algorithm: Beaumont et al. Biometrika 2009
#####################################################
.ABC_PMC<-function(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed=TRUE,seed_count=0,inside_prior=TRUE){
	T=length(tolerance_tab)
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)
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
	print("step 1 completed")
## steps 2 to T
	for (it in 2:T){
		nb_simul_step=nb_simul
		simul_below_tol2=NULL
		while (nb_simul_step>0){
			if (nb_simul_step>1){
				# Sampling of parameters around the previous particles
				tab_ini=.ABC_launcher_not_uniform_uni(model,prior_matrix,simul_below_tol[,1:nparam],tab_unfixed_param,tab_weight,nb_simul_step,use_seed,seed_count,inside_prior)
				seed_count=seed_count+nb_simul_step
				simul_below_tol2=rbind(simul_below_tol2,.selec_simul(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,tolerance_tab[it]))
				if (length(simul_below_tol2)>0){
					nb_simul_step=nb_simul-dim(simul_below_tol2)[1]
				}
			}
			else{
				tab_ini=.ABC_launcher_not_uniform_uni(model,prior_matrix,simul_below_tol[,1:nparam],tab_unfixed_param,tab_weight,nb_simul_step,use_seed,seed_count,inside_prior)
				seed_count=seed_count+nb_simul_step
				if (.compute_dist_single(summary_stat_target,tab_ini[(nparam+1):(nparam+nstat)],sd_simul)<tolerance_tab[it]){
					simul_below_tol2=rbind(simul_below_tol2,tab_ini)
					nb_simul_step=0
				}
			}
		} # until we get nb_simul simulations below the it-th tolerance threshold
		# update of particle weights
		tab_weight2=.compute_weight_uni(simul_below_tol2[,1:nparam][,tab_unfixed_param],simul_below_tol[,1:nparam][,tab_unfixed_param],tab_weight)
		# update of the set of particles and of the associated weights for the next ABC sequence
		tab_weight=tab_weight2
		simul_below_tol=matrix(0,nb_simul,(nparam+nstat))
		for (i1 in 1:nb_simul){
			for (i2 in 1:(nparam+nstat)){
				simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
			}
		}
		write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
		print(paste("step ",it," completed",sep=""))
	}
cbind(tab_weight,simul_below_tol)
}

## test
# linux
# .ABC_PMC(.binary_model("./parthy"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5),use_seed=TRUE,inside_prior=TRUE)
# windows
# .ABC_PMC(.binary_model("./parthy_test.exe"),prior_matrix,20,c(1,0.9,0.8),c(50,2.5),use_seed=TRUE,inside_prior=TRUE)


## function to select the alpha quantile closest simulations
############################################################
.selec_simul_alpha<-function(summary_stat_target,param,simul,sd_simul,alpha){
	dist=.compute_dist(summary_stat_target,simul,sd_simul)
	n_alpha=ceiling(alpha*length(dist))
	tol=sort(dist)[n_alpha]
	res=cbind(param[dist<=tol,],simul[dist<=tol,])
	if (dim(res)[1]==0){
		res=NULL
	}
res
}

## function to select the simulations that are at a distance smaller are equal to tol from the data
###################################################################################################
.selec_simulb<-function(summary_stat_target,param,simul,sd_simul,tol){
	dist=.compute_dist(summary_stat_target,simul,sd_simul)
	res=cbind(param[dist<=tol,],simul[dist<=tol,])
	if (dim(res)[1]==0){
		res=NULL
	}
res
}

## sequential algorithm of Drovandi & Pettitt 2011 - the proposal used is a multivariate normal (cf paragraph 2.2 - p. 227 in Drovandi & Pettitt 2011)
######################################################################################################################################################
.ABC_Drovandi<-function(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,alpha=0.5,c=0.01,first_tolerance_level_auto=TRUE,use_seed=TRUE,seed_count=0){
	n_alpha=ceiling(nb_simul*alpha)
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)
	tab_unfixed_param=array(TRUE,nparam)
	for (i in 1:nparam){
		tab_unfixed_param[i]=(prior_matrix[i,1]!=prior_matrix[i,2])
	}
	if (first_tolerance_level_auto){
		tol_end=tolerance_tab
	}
	else{
		tol_end=tolerance_tab[2]
	}

## step 1
	nb_simul_step=nb_simul
	simul_below_tol=NULL
	if (first_tolerance_level_auto){
		# classic ABC step
		tab_ini=ABC_rejection(model,prior_matrix,nb_simul_step,use_seed,seed_count)
		sd_simul=sd(tab_ini[,(nparam+1):(nparam+nstat)]) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
		seed_count=seed_count+nb_simul_step
		# selection of simulations below the first tolerance level
		simul_below_tol=rbind(simul_below_tol,.selec_simul_alpha(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,alpha))
		simul_below_tol=simul_below_tol[1:n_alpha,] # to be sure that there are not two or more simulations at a distance equal to the tolerance determined by the quantile
	}
	else{
	   nb_simul_step=n_alpha
	   while (nb_simul_step>0){
		if (nb_simul_step>1){
			# classic ABC step
			tab_ini=ABC_rejection(model,prior_matrix,nb_simul_step,use_seed,seed_count)
			sd_simul=sd(tab_ini[,(nparam+1):(nparam+nstat)]) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
			seed_count=seed_count+nb_simul_step
			# selection of simulations below the first tolerance level
			simul_below_tol=rbind(simul_below_tol,.selec_simul(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,tolerance_tab[1]))
			if (length(simul_below_tol)>0){
				nb_simul_step=n_alpha-dim(simul_below_tol)[1]
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
	  } # until we get n_alpha simulations below the first tolerance threshold
	}
	# initially, weights are equal
	tab_weight=array(1/n_alpha,n_alpha)
	write.table(cbind(tab_weight,simul_below_tol),file="output_step1",row.names=F,col.names=F,quote=F)
	print("step 1 completed")

## following steps until tol_end is reached
	tol_next=tolerance_tab[1]
	if (first_tolerance_level_auto){
		tol_next=max(.compute_dist(summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul))
	}
	R=1
	l=dim(simul_below_tol)[2]
	while (tol_next>tol_end){
		i_acc=0
		nb_simul_step=nb_simul-n_alpha
		simul_below_tol2=NULL
		for (i in 1:nb_simul_step){
			# pick a particle
			simul_picked=.particle_pick(simul_below_tol,tab_weight)
			for (j in 1:R){
				# move it
				param_moved=.move_particle(simul_picked[1:nparam][tab_unfixed_param],2*var(simul_below_tol[,1:nparam][,tab_unfixed_param]),prior_matrix[tab_unfixed_param,])
				param=simul_picked[1:nparam]
				param[tab_unfixed_param]=param_moved
				if (use_seed) {
					param=c((seed_count+i),param)
				}
				# perform a simulation
				new_simul=c(param,model(param))
				seed_count=seed_count+1
				if (use_seed) {
					new_simul=new_simul[2:(l+1)]
				}
				# check whether it is below tol_next and undo the move if it is not
				if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
					simul_picked=new_simul
					i_acc=i_acc+1
				}
			}
			simul_below_tol2=rbind(simul_below_tol2,simul_picked)
		}
		simul_below_tol2=rbind(simul_below_tol,simul_below_tol2)
		simul_below_tol=matrix(0,nb_simul,(nparam+nstat))
		for (i1 in 1:nb_simul){
			for (i2 in 1:(nparam+nstat)){
				simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
			}
		}
		tol_next=sort(.compute_dist(summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul))[n_alpha]
		simul_below_tol2=.selec_simulb(summary_stat_target,simul_below_tol[,1:nparam],simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul,tol_next) # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
		simul_below_tol=matrix(0,n_alpha,(nparam+nstat))
		for (i1 in 1:n_alpha){
			for (i2 in 1:(nparam+nstat)){
				simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
			}
		}
		p_acc=i_acc/(nb_simul_step*R)
		R=log(c)/log(1-p_acc)
	}

## final step to diversify the n_alpha particles
	simul_below_tol2=NULL
	for (i in 1:n_alpha){
		simul_picked=simul_below_tol[i,]
		for (j in 1:R){
			# move it
			param_moved=.move_particle(simul_picked[1:nparam][tab_unfixed_param],2*var(simul_below_tol[,1:nparam][,tab_unfixed_param]),prior_matrix[tab_unfixed_param,])
			param=simul_picked[1:nparam]
			param[tab_unfixed_param]=param_moved
			if (use_seed) {
				param=c((seed_count+i),param)
			}
			# perform a simulation
			new_simul=c(param,model(param))
			seed_count=seed_count+1
			if (use_seed) {
				new_simul=new_simul[2:(l+1)]
				}
			# check whether it is below tol_next and undo the move if it is not
			if (.compute_dist_single(summary_stat_target,as.numeric(new_simul[(nparam+1):(nparam+nstat)]),sd_simul)<=tol_next){ # we authorize the simulation to be equal to the tolerance level, for consistency with the quantile definition of the tolerance
				simul_picked=as.numeric(new_simul)
			}
		}
		simul_below_tol2=rbind(simul_below_tol2,simul_picked)
	}
cbind(tab_weight,simul_below_tol2)
}


## test
# linux
# .ABC_Drovandi(.binary_model("./parthy"),prior_matrix,20,c(1,0.8),c(50,2.5))
# .ABC_Drovandi(.binary_model("./parthy"),prior_matrix,20,c(1,0.8),c(50,2.5),first_tolerance_level_auto=FALSE)
# windows
# .ABC_Drovandi(.binary_model("./parthy_test.exe"),prior_matrix,20,c(1,0.8),c(50,2.5))
# .ABC_Drovandi(.binary_model("./parthy_test.exe"),prior_matrix,20,c(1,0.8),c(50,2.5),first_tolerance_level_auto=FALSE)

## rejection algorithm with M simulations per parameter set
############################################################
.ABC_rejection_M<-function(model,prior_matrix,nb_simul,M,use_seed,seed_count){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	for (i in 1:nb_simul){
		l=dim(prior_matrix)[1]
		param=array(0,l)
		for (j in 1:l){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		for (k in 1:M){
			if (use_seed) {
				param=c((seed_count+i),param)
				seed_count=seed_count+1
			}
			simul_summarystat=model(param)
			tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
			if (use_seed) {
				param=param[2:(l+1)]
			}
			tab_param=rbind(tab_param,param)
		}
	}
	cbind(tab_param,tab_simul_summarystat)
}

## function to compute a distance between a matrix of simulated statistics and the array of data summary statistics - for M replicates simulations
##################################################################################################################################################
.compute_dist_M<-function(M,summary_stat_target,simul,sd_simul){
	l=length(summary_stat_target)
	nsimul=dim(simul)[1]
	vartab=array(1,l)
	dist=array(0,nsimul)
	for (i in 1:l){
		vartab[i]=min(1,1/(sd_simul[i]*sd_simul[i])) ## differences between simul and data are normalized in each dimension by the empirical variances in each dimension
		dist=dist+vartab[i]*(simul[,i]-summary_stat_target[i])*(simul[,i]-summary_stat_target[i]) ## an euclidean distance is used
	}
	distb=matrix(0,nsimul/M,M)
	for (i in 1:(nsimul/M)){
		for (j in 1:M){
			distb[i,j]=dist[((i-1)*M+j)]
		}
	}
distb
}

## function to compute the updated weights, given a new tolerance value
#######################################################################
.compute_weight_delmoral<-function(particle_dist_mat,tolerance){
	n_particle=dim(particle_dist_mat)[1]
	new_weight=array(0,n_particle)
	for (i in 1:n_particle){
		new_weight[i]=length(particle_dist_mat[i,][particle_dist_mat[i,]<tolerance])
	}
new_weight/sum(new_weight)
}

## function to compute the updated ESS, given a new tolerance value
###################################################################
.compute_ESS<-function(particle_dist_mat,tolerance){
	n_particle=dim(particle_dist_mat)[1]
	new_weight=array(0,n_particle)
	for (i in 1:n_particle){
		new_weight[i]=length(particle_dist_mat[i,][particle_dist_mat[i,]<tolerance])
	}
	new_weight=new_weight/sum(new_weight)
1/(sum(new_weight*new_weight))
}

## function to compute the number of simul below a new tolerance value
######################################################################
.compute_below<-function(particle_dist_mat,tolerance){
	n_particle=dim(particle_dist_mat)[1]
	new_weight=array(0,n_particle)
	for (i in 1:n_particle){
		new_weight[i]=length(particle_dist_mat[i,][particle_dist_mat[i,]<tolerance])
	}
new_weight
}

## function to randomly pick a particle from a weighted array (of sum=1) for the Del Moral algorithm
####################################################################################################
.particle_pick_delmoral<-function(simul_below_tol,tab_weight,M){
	u=runif(1)
	weight_cum=cumsum(tab_weight)
	pos=1:length(tab_weight)
	p=min(pos[weight_cum>u])
	simul_below_tol[((1:M)+(p-1)*M),]
}

## function to replicate each cell of tab_weight M times
########################################################
.replicate_tab<-function(tab_weight,M){
	l=length(tab_weight)
	tab_weight2=array(0,M*l)
	for (i in 1:l){
		tab_weight2[((i-1)*M+(1:M))]=tab_weight[i]
	}
tab_weight2
}


## sequential algorithm of Del Moral et al. 2012 - the proposal used is a normal in each dimension (cf paragraph 3.2 in Del Moral et al. 2012)
##############################################################################################################################################
.ABC_Delmoral<-function(model,prior_matrix,nb_simul,alpha,M,nb_threshold,tolerance_target,summary_stat_target,use_seed=TRUE,seed_count=0){
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)
	tab_unfixed_param=array(TRUE,nparam)
	for (i in 1:nparam){
		tab_unfixed_param[i]=(prior_matrix[i,1]!=prior_matrix[i,2])
	}

# step 1
	# classic ABC step
	simul_below_tol=.ABC_rejection_M(model,prior_matrix,nb_simul,M,use_seed,seed_count)
	seed_count=seed_count+M*nb_simul
	tab_weight=rep(1/nb_simul,nb_simul)
	ESS=nb_simul
	uu=(1:nb_simul)*M  # to compute sd_simul with only one simulation per parameter set
	sd_simul=sd(simul_below_tol[uu,(nparam+1):(nparam+nstat)])  # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm
	l=dim(simul_below_tol)[2]
	if (M>1){
		particle_dist_mat=.compute_dist_M(M,summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul)
	}
	else{
		particle_dist_mat=.compute_dist(summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul)
	}
	dim(particle_dist_mat)<-c(nb_simul,M)
	new_tolerance=max(particle_dist_mat)

	tab_weight2=.replicate_tab(tab_weight,M)
	write.table(cbind(tab_weight2,simul_below_tol),file="output_step1",row.names=F,col.names=F,quote=F)
	print("step 1 completed")

# following steps
	kstep=1
   while(new_tolerance>tolerance_target){	
	kstep=kstep+1
	# determination of the new tolerance
	ESS_target=alpha*ESS

	tolerance_list=sort(as.numeric(names(table(particle_dist_mat))),dec=TRUE)
	i=1
	test=FALSE
	while((!test)&&(i<length(tolerance_list))){
		i=i+1
		# computation of new ESS with the new tolerance value
		new_ESS=.compute_ESS(particle_dist_mat,tolerance_list[i])
		# check whether this value is below ESS_targ
		if (new_ESS<ESS_target){
			new_tolerance=tolerance_list[(i-1)]
			test=TRUE
		}
	}

	# if effective sample size is too small, resampling of particles
	ESS=.compute_ESS(particle_dist_mat,new_tolerance)
	tab_weight=.compute_weight_delmoral(particle_dist_mat,new_tolerance)
	tab_below=.compute_below(particle_dist_mat,new_tolerance)
	particles=matrix(0,(nb_simul*M),(nparam+nstat))
	if (ESS<nb_threshold){
		# sample nb_simul particles 
		for (i in 1:nb_simul){
			particles[((1:M)+(i-1)*M),]=as.matrix(.particle_pick_delmoral(simul_below_tol,tab_weight,M))
		}
		simul_below_tol=matrix(0,nb_simul*M,(nparam+nstat))
		for (i1 in 1:(nb_simul*M)){
			for (i2 in 1:(nparam+nstat)){
				simul_below_tol[i1,i2]=as.numeric(particles[i1,i2])
			}
		}
		particles=particles[uu,1:nparam]
		if (M>1){
			particle_dist_mat=.compute_dist_M(M,summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul)
		}
		else{
			particle_dist_mat=.compute_dist(summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul)
		}
		dim(particle_dist_mat)<-c(nb_simul,M)
		tab_below=.compute_below(particle_dist_mat,new_tolerance)
		# reset their weight to 1/nb_simul
		tab_weight=rep(1/nb_simul,nb_simul)
		ESS=nb_simul
	}
	else{
		particles=simul_below_tol[uu,1:nparam]
	}

	# MCMC move
	covmat=2*cov.wt(particles[,tab_unfixed_param][tab_weight>0,],as.vector(tab_weight[tab_weight>0]))$cov
	l_array=dim(particles[,tab_unfixed_param])[2]
	sd_array=array(1,l_array)
	for (j in 1:l_array){
		sd_array[j]=sqrt(covmat[j,j])

	}
	simul_below_tol2=simul_below_tol
	simul_below_tol=matrix(0,nb_simul*M,(nparam+nstat))
	for (i in 1:nb_simul){
		if (tab_weight[i]>0){
			tab_new_simul=NULL
			# move it
			param_moved=.move_particle_uni(as.numeric(particles[i,tab_unfixed_param]),sd_array,prior_matrix[tab_unfixed_param,])
			param=particles[i,]
			param[tab_unfixed_param]=param_moved
			if (use_seed) {
				param=c((seed_count+i),param)
			}
			# perform M simulations
			for (j in 1:M){
				new_simul=c(param,model(param))
				seed_count=seed_count+1
				if (use_seed) {
					new_simul=new_simul[2:(l+1)]
				}
				tab_new_simul=rbind(tab_new_simul,new_simul)
			}
			if (M>1){
				tab_new_simul2=matrix(0,M,(nparam+nstat))
				for (i1 in 1:M){
					for (i2 in 1:(nparam+nstat)){
						tab_new_simul2[i1,i2]=as.numeric(tab_new_simul[i1,i2])
					}
				}
			}
			else{
				tab_new_simul2=as.numeric(tab_new_simul)
			}
			dim(tab_new_simul2)<-c(M,(nparam+nstat))
			# check whether the move is accepted
			n_acc=1
			if (M>1){
				new_dist=.compute_dist_M(M,summary_stat_target,tab_new_simul2[,(nparam+1):(nparam+nstat)],sd_simul)
				n_acc=length(new_dist[new_dist<new_tolerance])
			}
			else{
				new_dist=.compute_dist(summary_stat_target,rbind(tab_new_simul2[(nparam+1):(nparam+nstat)],tab_new_simul2[(nparam+1):(nparam+nstat)]),sd_simul)
				if (new_dist[1]>new_tolerance){
					n_acc=0
				}
			}
			MH=min(1,(n_acc/tab_below[i]))
			uuu=runif(1)
			if (uuu<=MH){
				for (i1 in 1:M){
					for (i2 in 1:(nparam+nstat)){
						simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(tab_new_simul2[i1,i2])
					}
				}

			}
			else{
				for (i1 in 1:M){
					for (i2 in 1:(nparam+nstat)){
						simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(simul_below_tol2[((i1)+(i-1)*M),i2])
					}
				}
			}
		}
		else{
			for (i1 in 1:M){
				for (i2 in 1:(nparam+nstat)){
					simul_below_tol[((i1)+(i-1)*M),i2]=as.numeric(simul_below_tol2[((i1)+(i-1)*M),i2])
				}
			}
		}
	}	
	if (M>1){
		particle_dist_mat=.compute_dist_M(M,summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul)
	}
	else{
		particle_dist_mat=.compute_dist(summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul)
	}
	dim(particle_dist_mat)<-c(nb_simul,M)
	tab_weight=.compute_weight_delmoral(particle_dist_mat,new_tolerance)
	tab_weight2=.replicate_tab(tab_weight,M)
	write.table(cbind(tab_weight2,simul_below_tol),file=paste("output_step",kstep,sep=""),row.names=F,col.names=F,quote=F)
	print(paste("step ",kstep," completed",sep=""))
   }
	
cbind(tab_weight2,simul_below_tol)
}


## test
# linux
# .ABC_Delmoral(.binary_model("./parthy"),prior_matrix,20,0.5,1,5,0.8,c(50,2.5))
# .ABC_Delmoral(.binary_model("./parthy"),prior_matrix,20,0.5,4,5,0.8,c(50,2.5))
# windows
# .ABC_Delmoral(.binary_model("./parthy_test.exe"),prior_matrix,10,0.5,1,3,0.8,c(50,2.5))
# .ABC_Delmoral(.binary_model("./parthy_test.exe"),prior_matrix,20,0.5,4,5,0.8,c(50,2.5))


library(lhs)

## function to sample in the prior distributions using a Latin Hypercube sample
###############################################################################
.ABC_rejection_lhs<-function(model,prior_matrix,nb_simul,tab_unfixed_param,use_seed=TRUE,seed_count=0){
	tab_simul_summarystat=NULL
	tab_param=NULL
	l=dim(prior_matrix)[1]
	nparam=length(tab_unfixed_param[tab_unfixed_param])
	random_tab=randomLHS(nb_simul,nparam)

	for (i in 1:nb_simul){
		param=prior_matrix[,1]
		for (j in 1:nparam){
			param[tab_unfixed_param][j]=prior_matrix[tab_unfixed_param,][j,1]+(prior_matrix[tab_unfixed_param,][j,2]-prior_matrix[tab_unfixed_param,][j,1])*random_tab[i,j]
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

## function to compute particle weights without normalizing to 1
################################################################
.compute_weightb<-function(param_simulated,param_previous_step,tab_weight,prior_density){
	vmat=2*var(param_previous_step)
	n_particle=dim(param_previous_step)[1]
	n_new_particle=dim(param_simulated)[1]
	tab_weight_new=array(0,n_new_particle)
	for (i in 1:n_particle){
		for (j in 1:n_new_particle){
			tab_weight_new[j]=tab_weight_new[j]+tab_weight[i]*dmnorm(param_simulated[j,],param_previous_step[i,],vmat)
		}
	}
	tab_weight_new=prior_density/tab_weight_new
tab_weight_new
}

## function to perform ABC simulations from a non-uniform prior (derived from a set of particles)
#################################################################################################
.ABC_launcher_not_uniformc<-function(model,prior_matrix,param_previous_step,tab_unfixed_param,tab_weight,nb_simul,use_seed,seed_count,inside_prior){
	tab_simul_summarystat=NULL
	tab_param=NULL
	k_acc=0
	for (i in 1:nb_simul){
		l=dim(param_previous_step)[2]
		if (!inside_prior){
			k_acc=k_acc+1
			# pick a particle
			param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
			# move it
			param_moved=.move_particle(param_picked,2*cov.wt(param_previous_step[,tab_unfixed_param],as.vector(tab_weight))$cov,prior_matrix[tab_unfixed_param,]) # only variable parameters are moved, computation of a WEIGHTED variance
		}
		else{
			test=FALSE
			while (!test){
				k_acc=k_acc+1
				# pick a particle
				param_picked=.particle_pick(param_previous_step[,tab_unfixed_param],tab_weight)
				# move it
				param_moved=.move_particle(param_picked,2*cov.wt(param_previous_step[,tab_unfixed_param],as.vector(tab_weight))$cov,prior_matrix[tab_unfixed_param,]) # only variable parameters are moved, computation of a WEIGHTED variance
				test=.is_included(param_moved,prior_matrix)
			}

		}
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
list(cbind(tab_param,tab_simul_summarystat),nb_simul/k_acc)
}

## sequential algorithm of Lenormand et al. 2012 
################################################
.ABC_Lenormand<-function(model,prior_matrix,nb_simul,summary_stat_target,alpha=0.5,p_acc_min=0.05,use_seed=TRUE,seed_count=0,inside_prior=TRUE){
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)
	tab_unfixed_param=array(TRUE,nparam)
	for (i in 1:nparam){
		tab_unfixed_param[i]=(prior_matrix[i,1]!=prior_matrix[i,2])
	}
	n_alpha=ceiling(nb_simul*alpha)
	prior_density=1
	for (i in 1:nparam){
		if (tab_unfixed_param[i]){
			prior_density=prior_density/(prior_matrix[i,2]-prior_matrix[i,1])
		}
	}

## step 1
	# ABC rejection step with LHS
	tab_ini=.ABC_rejection_lhs(model,prior_matrix,nb_simul,tab_unfixed_param,use_seed,seed_count)
	seed_count=seed_count+nb_simul
	sd_simul=sd(tab_ini[,(nparam+1):(nparam+nstat)]) # determination of the normalization constants in each dimension associated to each summary statistic, this normalization will not change during all the algorithm

	# selection of the alpha quantile closest simulations
	simul_below_tol=NULL
	simul_below_tol=rbind(simul_below_tol,.selec_simul_alpha(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],sd_simul,alpha))
	simul_below_tol=simul_below_tol[1:n_alpha,] # to be sure that there are not two or more simulations at a distance equal to the tolerance determined by the quantile

	# initially, weights are equal
	tab_weight=array(1,n_alpha)

	write.table(cbind(tab_weight,simul_below_tol),file="output_step1",row.names=F,col.names=F,quote=F)
	print("step 1 completed")
	tab_dist=.compute_dist(summary_stat_target,simul_below_tol[,(nparam+1):(nparam+nstat)],sd_simul)
	tol_next=max(tab_dist)

## following steps
	p_acc=p_acc_min+1
	nb_simul_step=nb_simul-n_alpha
	it=1
	while (p_acc>p_acc_min){
		it=it+1
		simul_below_tol2=NULL
		tab_inic=.ABC_launcher_not_uniformc(model,prior_matrix,simul_below_tol[,1:nparam],tab_unfixed_param,tab_weight/sum(tab_weight),nb_simul_step,use_seed,seed_count,inside_prior)
		tab_ini=tab_inic[[1]]
		seed_count=seed_count+nb_simul_step
		if (!inside_prior){
			tab_weight2=.compute_weightb(tab_ini[,1:nparam][,tab_unfixed_param],simul_below_tol[,1:nparam][,tab_unfixed_param],tab_weight/sum(tab_weight),prior_density)
		}
		else{
			tab_weight2=tab_inic[[2]]*(.compute_weightb(tab_ini[,1:nparam][,tab_unfixed_param],simul_below_tol[,1:nparam][,tab_unfixed_param],tab_weight/sum(tab_weight),prior_density))
		}
		simul_below_tol2=rbind(simul_below_tol,as.matrix(tab_ini))
		tab_weight=c(tab_weight,tab_weight2)
		tab_dist2=.compute_dist(summary_stat_target,tab_ini[,(nparam+1):(nparam+nstat)],sd_simul)
		p_acc=length(tab_dist2[tab_dist2<=tol_next])/nb_simul_step
		tab_dist=c(tab_dist,tab_dist2)
		tol_next=sort(tab_dist)[n_alpha]
		simul_below_tol2=simul_below_tol2[tab_dist<=tol_next,]
		tab_weight=tab_weight[tab_dist<=tol_next]
		tab_weight=tab_weight[1:n_alpha]
		tab_dist=tab_dist[tab_dist<=tol_next]
		tab_dist=tab_dist[1:n_alpha]
		simul_below_tol=matrix(0,n_alpha,(nparam+nstat))
		for (i1 in 1:n_alpha){
			for (i2 in 1:(nparam+nstat)){
				simul_below_tol[i1,i2]=as.numeric(simul_below_tol2[i1,i2])
			}
		}
		write.table(as.matrix(cbind(tab_weight,simul_below_tol)),file=paste("output_step",it,sep=""),row.names=F,col.names=F,quote=F)
		print(paste("step ",it," completed",sep=""))
	}
cbind(tab_weight,simul_below_tol)
}

## test
# linux
# .ABC_Lenormand(.binary_model("./parthy"),prior_matrix,20,c(50,2.5),inside_prior=TRUE)
# windows
# .ABC_Lenormand(.binary_model("./parthy_test.exe"),prior_matrix,20,c(50,2.5),inside_prior=TRUE)


## function to move a particle with a unidimensional uniform jump
#################################################################
.move_particle_uni_uniform<-function(param_picked,sd_array,prior_matrix){
	test=FALSE
	res=param_picked
	while (!test){
		for (i in 1:length(param_picked)){
			res[i]=runif(n = 1, min = param_picked[i]-sd_array[i],max=param_picked[i]+sd_array[i])
		}
		test=.is_included(res,prior_matrix)
	}
res
}

## ABC-MCMC algorithm of Marjoram et al. 2003
#############################################
.ABC_MCMC<-function(model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,dist_max,tab_normalization,proposal_range,use_seed=TRUE,seed_count=0){
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)
	tab_simul_summary_stat=NULL
	tab_param=NULL
	tab_unfixed_param=array(TRUE,nparam)
	for (i in 1:nparam){
		tab_unfixed_param[i]=(prior_matrix[i,1]!=prior_matrix[i,2])
	}

	# initial draw of a particle below the tolerance dist_max
	test=FALSE
	while (!test){
		param=array(0,nparam)
		for (j in 1:nparam){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		if (use_seed) {
			param=c((seed_count+i),param)
		}
		simul_summary_stat=model(param)
		dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),tab_normalization)
		if (dist_simul<dist_max){
			test=TRUE
		}
		seed_count=seed_count+1
	}
	tab_simul_summary_stat=rbind(tab_simul_summary_stat,simul_summary_stat)
	tab_param=rbind(tab_param,param)
	if (use_seed) {
			tab_param=tab_param[,2:(nparam+1)]
	}
	tab_simul_ini=as.numeric(simul_summary_stat)
	param_ini=tab_param
	print("initial draw performed ")

	# chain run
	tab_param=param_ini
	tab_simul_summary_stat=tab_simul_ini
	for (is in 2:n_obs){
		for (i in 1:n_between_sampling){
			param=.move_particle_uni_uniform(as.numeric(param_ini),proposal_range,prior_matrix)
			if (use_seed) {
				param=c(seed_count,param)
			}	
			simul_summary_stat=model(param)
			if (use_seed) {
				param=param[2:(nparam+1)]
			}
			dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),tab_normalization)
			if (dist_simul<dist_max){
				param_ini=param
				tab_simul_ini=as.numeric(simul_summary_stat)
			}
			seed_count=seed_count+1
		}
		tab_simul_summary_stat=rbind(tab_simul_summary_stat,tab_simul_ini)
		tab_param=rbind(tab_param,as.numeric(param_ini))
		if (is%%100==0){
			print(paste(is," ",sep=""))
		}
	}	
cbind(tab_param,tab_simul_summary_stat)	
}

## test
# linux
# .ABC_MCMC(.binary_model("./parthy"),prior_matrix,22,10,c(50,2.5),1,c(33,1),c(25,25,1,0,0))
# windows
# .ABC_MCMC(.binary_model("./parthy_test.exe"),prior_matrix,22,10,c(50,2.5),1,c(33,1),c(25,25,1,0,0))

## ABC-MCMC algorithm of Marjoram et al. 2003 with automatic determination of the tolerance and proposal range following Wegmann et al. 2009
############################################################################################################################################
.ABC_MCMC2<-function(model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,n_calibration=10000,tolerance_quantile=0.01,proposal_phi=1,use_seed=TRUE,seed_count=0){
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)
	tab_simul_summary_stat=NULL
	tab_param=NULL
	tab_unfixed_param=array(TRUE,nparam)
	for (i in 1:nparam){
		tab_unfixed_param[i]=(prior_matrix[i,1]!=prior_matrix[i,2])
	}

	# initial draw of a particle
	for (i in 1:(n_calibration)){
		param=array(0,nparam)
		for (j in 1:nparam){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		if (use_seed) {
			param=c((seed_count+i),param)
		}
		simul_summary_stat=model(param)
		tab_simul_summary_stat=rbind(tab_simul_summary_stat,simul_summary_stat)
		tab_param=rbind(tab_param,param)
	}
	if (use_seed) {
			tab_param=tab_param[,2:(nparam+1)]
	}
	proposal_range=array(0,nparam)
	for (i in 1:nparam){
		proposal_range[i]=sd(tab_param[,i])*proposal_phi
	}
	simuldist=.compute_dist(summary_stat_target,tab_simul_summary_stat,proposal_range/proposal_phi)
	ord_sim=order(simuldist,decreasing=F)
	nmax=ceiling(tolerance_quantile*n_calibration)
	dist_max=simuldist[(ord_sim[nmax])]
	tab_param=tab_param[(ord_sim[1:nmax]),]
	n_ini=sample(nmax,1)
	tab_simul_ini=as.numeric(tab_simul_summary_stat[(ord_sim[n_ini]),])
	param_ini=tab_param[n_ini,]
	print("initial calibration performed ")

	# chain run
	tab_param=param_ini
	tab_simul_summary_stat=tab_simul_ini
	seed_count=seed_count+n_calibration+1
	for (is in 2:n_obs){
		for (i in 1:n_between_sampling){
			param=.move_particle_uni_uniform(as.numeric(param_ini),proposal_range,prior_matrix)
			if (use_seed) {
				param=c(seed_count,param)
			}	
			simul_summary_stat=model(param)
			if (use_seed) {
				param=param[2:(nparam+1)]
			}
			dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),proposal_range/proposal_phi)
			if (dist_simul<dist_max){
				param_ini=param
				tab_simul_ini=as.numeric(simul_summary_stat)
			}
			seed_count=seed_count+1
		}
		tab_simul_summary_stat=rbind(tab_simul_summary_stat,tab_simul_ini)
		tab_param=rbind(tab_param,as.numeric(param_ini))
		if (is%%100==0){
			print(paste(is," ",sep=""))
		}
	}	
cbind(tab_param,tab_simul_summary_stat)	
}

## test
# linux
# .ABC_MCMC2(.binary_model("./parthy"),prior_matrix,22,10,c(50,2.5),n_calibration=10,tolerance_quantile=0.2,proposal_phi=1)
# windows
# .ABC_MCMC2(.binary_model("./parthy_test.exe"),prior_matrix,22,10,c(50,2.5),n_calibration=10,tolerance_quantile=0.2,proposal_phi=1)

library(pls)
library(MASS)

## ABC-MCMC algorithm of Wegmann et al. 2009 - the PLS step is drawn from the manual of ABCtoolbox (figure 9) - NB: for consistency with ABCtoolbox, AM11-12 are not implemented in the algorithm
#################################################################################################################################################################################################
.ABC_MCMC3<-function(model,prior_matrix,n_obs,summary_stat_targ,n_calibration=10000,tolerance_quantile=0.01,proposal_phi=1,n_between_sampling=1,numcomp=0,use_seed=TRUE,seed_count=0){
##AM1
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_targ)
	if (numcomp==0){
		numcomp=nstat
	}
	tab_simul_summary_stat=NULL
	tab_param=NULL
	tab_unfixed_param=array(TRUE,nparam)
	for (i in 1:nparam){
		tab_unfixed_param[i]=(prior_matrix[i,1]!=prior_matrix[i,2])
	}

	# initial draw of a particle
	for (i in 1:(n_calibration)){
		param=array(0,nparam)
		for (j in 1:nparam){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		if (use_seed) {
			param=c((seed_count+i),param)
		}
		simul_summary_stat=model(param)
		tab_simul_summary_stat=rbind(tab_simul_summary_stat,simul_summary_stat)
		tab_param=rbind(tab_param,param)
	}
	if (use_seed) {
			tab_param=tab_param[,2:(nparam+1)]
	}

## AM2: PLS step
	print("AM2 ")
	#standardize the params
	sparam=tab_param[,tab_unfixed_param]
	ls=dim(sparam)[2]
	for(i in 1:ls){
		sparam[,i]=(sparam[,i]-mean(sparam[,i]))/sd(sparam[,i])
	}
	#force stat in [1,2]
	myMax<-c()
	myMin<-c()
	lambda<-c()
	myGM<-c()
	stat=tab_simul_summary_stat
	print("stat 1 ")
	print(stat)
	summary_stat_target=summary_stat_targ
	for (i in 1:nstat){
		myMax<-c(myMax,max(stat[,i]))
		myMin<-c(myMin,min(stat[,i]))
		stat[,i]=1+(stat[,i]-myMin[i])/(myMax[i]-myMin[i])
		summary_stat_target[i]=1+(summary_stat_target[i]-myMin[i])/(myMax[i]-myMin[i])
	}
	print("stat 2 ")
	print(stat)
	#transform statistics via boxcox
	dmat=matrix(0,n_calibration,(ls+1))
	for(i in 1:nstat){
		d=cbind(as.vector(as.numeric(stat[,i])),as.matrix(sparam))
		for (i1 in 1:n_calibration){
			for (i2 in 1:(ls+1)){
				dmat[i1,i2]=as.numeric(d[i1,i2])
			}
		}
		mylm<-lm(as.formula(as.data.frame(dmat)),data=as.data.frame(dmat))
		#mylm<-lm(stat[,i]~as.matrix(sparam))
		myboxcox<-boxcox(mylm,lambda=seq(-20,100,1/10),interp=T,eps=1/50)
		lambda<-c(lambda,myboxcox$x[myboxcox$y==max(myboxcox$y)])
		myGM<-c(myGM, exp(mean(log(stat[,i]))))
	}
	#standardize the BC-stat
	myBCMeans<-c()
	myBCSDs<-c()
	for(i in 1:nstat){
		stat[,i]<-((stat[,i]^lambda[i]) - 1)/(lambda[i]*(myGM[i]^(lambda[i]-1)));
		summary_stat_target[i]<-((summary_stat_target[i]^lambda[i]) - 1)/(lambda[i]*(myGM[i]^(lambda[i]-1)));
		myBCSDs<-c(myBCSDs, sd(stat[,i]))
		myBCMeans<-c(myBCMeans, mean(stat[,i]))
		stat[,i]<-(stat[,i]-myBCMeans[i])/myBCSDs[i]
		summary_stat_target[i]<-(summary_stat_target[i]-myBCMeans[i])/myBCSDs[i]
	}
	#perform pls
	myPlsr<-plsr(as.matrix(sparam)~as.matrix(stat), scale=F,ncomp=numcomp,validation='LOO')
	pls_transformation=matrix(0,numcomp,nstat)
	for (i in 1:numcomp){
		pls_transformation[i,]=as.numeric(myPlsr$loadings[,i])
	}

## AM3
	print("AM3 ")
	summary_stat_target=t(pls_transformation %*% t(summary_stat_target))
	stat_pls=t(pls_transformation %*% t(stat))
	simuldist=.compute_dist(summary_stat_target,stat_pls,rep(1,numcomp))

## AM4
	print("AM4 ")
	ord_sim=order(simuldist,decreasing=F)
	nmax=ceiling(tolerance_quantile*n_calibration)
	dist_max=simuldist[(ord_sim[nmax])]

	proposal_range=array(0,nparam)
	for (i in 1:nparam){
		proposal_range[i]=sd(tab_param[,i])*proposal_phi
	}
	tab_param=tab_param[(ord_sim[1:nmax]),]
	print("initial calibration performed ")

## AM5: chain run
	print("AM5 ")
	n_ini=sample(nmax,1)
	tab_simul_ini=as.numeric(tab_simul_summary_stat[(ord_sim[n_ini]),])
	param_ini=tab_param[n_ini,]
	tab_param=param_ini
	tab_simul_summary_stat=tab_simul_ini
	seed_count=seed_count+n_calibration+1
	for (is in 2:n_obs){
		for (i in 1:n_between_sampling){
## AM6
	print("AM6 ")
			param=.move_particle_uni_uniform(as.numeric(param_ini),proposal_range,prior_matrix)
			if (use_seed) {
				param=c(seed_count,param)
			}
## AM7	
	print("AM7 ")
			simul_summary_stat=model(param)
			if (use_seed) {
				param=param[2:(nparam+1)]
			}
			for (ii in 1:nstat){
				simul_summary_stat[ii]=1+(simul_summary_stat[ii]-myMin[ii])/(myMax[ii]-myMin[ii])
			}
			for(ii in 1:nstat){
				simul_summary_stat[ii]<-(simul_summary_stat[ii]-myBCMeans[ii])/myBCSDs[ii]
			}
			simul_summary_stat=t(pls_transformation %*% t(simul_summary_stat))
			dist_simul=.compute_dist_single(summary_stat_target,as.numeric(simul_summary_stat),rep(1,numcomp))
## AM8-9
	print("AM8-9 ")
			if (dist_simul<dist_max){
				param_ini=param
				tab_simul_ini=as.numeric(simul_summary_stat)
			}
			seed_count=seed_count+1
		}
		tab_simul_summary_stat=rbind(tab_simul_summary_stat,tab_simul_ini)
		tab_param=rbind(tab_param,as.numeric(param_ini))
		if (is%%100==0){
			print(paste(is," ",sep=""))
		}
	}	
cbind(tab_param,tab_simul_summary_stat)	
}

## test
# linux
# .ABC_MCMC3(.binary_model("./parthy"),prior_matrix,22,c(50,2.5),n_calibration=10,tolerance_quantile=0.2,proposal_phi=1)
# windows
# .ABC_MCMC3(.binary_model("./parthy_test.exe"),prior_matrix,22,c(50,2.5),n_calibration=10,tolerance_quantile=0.2,proposal_phi=1)
