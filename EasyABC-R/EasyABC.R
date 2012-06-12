ABC_launcher<-function(model,prior_matrix,nb_simul){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	
	for (i in 1:nb_simul){
		l=dim(prior_matrix)[1]
		param=array(0,l)
		for (j in 1:l){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		write.table(param,file="input",row.names=F,col.names=F,quote=F)
		system(model)
		simul_summarystat=read.table("output",h=F)
		tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
		tab_param=rbind(tab_param,param)
	}
cbind(tab_param,tab_simul_summarystat)
}


ABC_launcher_with_seed<-function(model,prior_matrix,nb_simul){
	tab_simul_summarystat=NULL
	tab_param=NULL
	
	
	for (i in 1:nb_simul){
		l=dim(prior_matrix)[1]
		param=array(0,l)
		for (j in 1:l){
			param[j]=runif(1,min=prior_matrix[j,1],max=prior_matrix[j,2])
		}
		param=c(i,param)
		write.table(param,file="input",row.names=F,col.names=F,quote=F)
		system(model)
		simul_summarystat=read.table("output",h=F)
		tab_simul_summarystat=rbind(tab_simul_summarystat,simul_summarystat)
		tab_param=rbind(tab_param,param)
	}
cbind(tab_param,tab_simul_summarystat)
}



prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
prior_matrix

ABC_launcher_with_seed("parthy_test.exe",prior_matrix,10)
