ABC_launcher<-function(model,prior_matrix,nb_simul,use_seed=FALSE){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	for (i in 1:nb_simul){
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
	cbind(tab_param,tab_simul_summarystat)
}

binary_model<-function(command) {
	invoke<-function(param) {
		write.table(param,file="input",row.names=F,col.names=F,quote=F)
		system(command)
		read.table("output",h=F)
	}
}

# sample usage

prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
prior_matrix

ABC_launcher(binary_model("./parthy"),prior_matrix,10,TRUE)


## PMC ABC algorithm: Beaumont et al. Biometrika 2009

compute_dist<-function(summary_stat_target,simul){
	l=length(summary_stat_target)
	nsimul=dim(simul)[1]
	vartab=array(1,l)
	dist=array(0,nsimul)
	for (i in 1:l){
		vartab[i]=min(1,1/var(simul[,i]))
		dist=dist+vartab[i]*(simul[,i]-summary_stat_target[i])*(simul[,i]-summary_stat_target[i])
	}
dist
}

compute_dist_single<-function(summary_stat_target,simul,sd_simul){
	l=length(summary_stat_target)
	dist=0
	for (i in 1:l){
		vartab[i]=min(1,1/(sd_simul[i]*sd_simul[i]))
		dist=dist+vartab[i]*(simul[i]-summary_stat_target[i])*(simul[i]-summary_stat_target[i])
	}
dist
}

selec_simul<-function(summary_stat_target,param,simul,tol){
	dist=compute_dist(summary_stat_target,simul)
cbind(param[dist<tol],simul[dist<tol])
}

particle_pick<-function(tab_weight){
	u=runif(1)
	weight_cum=cumsum(tab_weight)
	pos=1:length(tab_weight)
	p=min(pos[weight_cum>u])
param[p,]
}

# with package mnormt
library(mnormt)

move_particle<-function(param_picked,varcov_matrix){
rmnorm(n = 1, mean = param_picked, varcov_matrix)
}

ABC_launcher_not_uniform<-function(model,param_previous_step,tab_weight,nb_simul,use_seed){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	for (i in 1:nb_simul){
		l=dim(prior_matrix)[1]
		# pick a particle
		param_picked=particle_pick(tab_weight)
		# move it
		param=move_particle(param_picked,2*var(param_previous_step))
		if (use_seed) {
			param=c(i,param)
		}
		simul_summarystat=model(param)
		tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
		tab_param=rbind(tab_param,param)
	}
	cbind(tab_param,tab_simul_summarystat)
}

compute_weight<-function(param_simulated,param_previous_step,tab_weight){
	vmat=2*var(param_previous_step)
	n_particle=dim(param_previous_step)[1]
	n_new_particle=dim(param_simulated_step)[1]
	tab_weight_new=array(0,n_new_particle)
	for (i in 1:n_particle){
		for (j in 1:n_new_particle){
			tab_weight_new[j]=tab_weight_new[j]+tab_weight[i]*dmnorm(param_simulated[j,],param_previous_step[i,],vmat)
		}
	}
tab_weight_new/sum(tab_weight_new)
}

ABC_PMC<-function(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed=FALSE){
	T=length(epsilon_tab)
	nparam=dim(prior_matrix)[1]
	nstat=length(summary_stat_target)

## step 1
	nb_simul_step=nb_simul
	simul_below_tol=NULL
	while (nb_simul_step>0){
		tab_ini=ABC_launcher(model,prior_matrix,nb_simul_step,use_seed)
		if (use_seed){
			tab_ini=tab_ini[,2:(dim(tab_ini)[2])]
		}
		
		simul_below_tol=rbind(simul_below_tol,selec_simul(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],tolerance_tab[1]))
		if (length(simul_below_tol)>0){
			nb_simul_step=nb_simul-dim(simul_below_tol)[1]
		}
	}
	tab_weight=array(1/nb_simul,nb_simul)

## steps 2 to T
	for (it in 2:T){
		simul_below_tol2=NULL
		while (nb_simul_step>0){
			tab_ini=ABC_launcher_not_uniform(model,simul_below_tol[,1:nparam],tab_weight,nb_simul_step,use_seed)
			if (use_seed){
				tab_ini=tab_ini[,2:(dim(tab_ini)[2])]
			}
			simul_below_tol2=rbind(simul_below_tol2,selec_simul(summary_stat_target,tab_ini[,1:nparam],tab_ini[,(nparam+1):(nparam+nstat)],tolerance_tab[it]))
			if (length(simul_below_tol2)>0){
				nb_simul_step=nb_simul-dim(simul_below_tol2)[1]
			}
		}
		tab_weight2=compute_weight(simul_below_tol2[,1:nparam],simul_below_tol[,1:nparam],tab_weight)
		tab_weight=tab_weight2
		simul_below_tol=simul_below_tol2
	}
cbind(tab_weight,simul_below_tol)
}

## Algo de Marjoram

## Algo de Marjoram avec détermination de epsilon et proposal_range de Wegmann et al. 2009
ABC_MCMC2<-function(model,prior_matrix,n_obs,n_between_sampling,summary_stat_target,proposal_method="unif",burn_in_length=0,n_calibration=10000,tolerance_quantile=0.01,proposal_phi=1,use_seed=FALSE){
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
	simuldist=compute_dist(summary_stat_target,tab_simul_summarystat)
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
			param=move_particle(param_ini,proposal_method,proposal_range,use_seed,seed_compt)
			simul_summarystat=model(param)
			dist_simul=compute_dist_single(summary_stat_target,simul_summary_stat,proposal_range/phi)
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

