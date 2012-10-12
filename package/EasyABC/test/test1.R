# library nécessaire : mnormt


sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
# pour eviter les effets de bords
rm(list = ls())
sourceDir('R')



# brute-force ABC
# .ABC_PMC2(.binary_model("./parthy"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5),use_seed=TRUE,inside_prior=TRUE)
prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
prior_matrix
tmps.ini=Sys.time()         #top chrono
ABC_rejection(.binary_model("../../model/parthy/parthy"),prior_matrix,10,TRUE)
cat("temps mis : ",difftime(Sys.time(),tmps.ini,units="sec")," s\n") #top fin chrono


# sequential algorithm
# .ABC_PMC(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed,seed_count,inside_prior),
.ABC_PMC(.binary_model("../../model/parthy/parthy"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5),use_seed=TRUE,inside_prior=TRUE)

#PMC2 = .ABC_PMC2(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,use_seed,seed_count,inside_prior),
ABC_sequential('PMC2',.binary_model("../../model/parthy/parthy"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5))

#Drovandi = .ABC_Drovandi(model,prior_matrix,nb_simul,tolerance_tab,summary_stat_target,alpha,c,first_tolerance_level_auto,use_seed,seed_count),
ABC_sequential('Drovandi',.binary_model("../../model/parthy/parthy"),prior_matrix,20,c(0.8,0.6,0.4),c(50,2.5))

#Delmoral = .ABC_Delmoral(model,prior_matrix,nb_simul,alpha,M,nb_threshold,tolerance_target,summary_stat_target,use_seed,seed_count),

#Lenormand = .ABC_Lenormand
