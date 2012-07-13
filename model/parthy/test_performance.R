source("../../EasyABC-R/EasyABC.R")
prior_matrix=c(1,1,-1,1000,10000,100,100,1,1000,10000)
dim(prior_matrix)<-c(5,2)
prior_matrix
# linux
tmps.ini=Sys.time()         #top chrono
ABC_rejection(.binary_model("./parthy"),prior_matrix,10,TRUE)
cat("temps mis : ",difftime(Sys.time(),tmps.ini,units="sec")," s\n") #top fin chrono

tmps.ini=Sys.time()         #top chrono
ABC_rejection(.binary_model("./parthy"),prior_matrix,100,TRUE)
cat("temps mis : ",difftime(Sys.time(),tmps.ini,units="sec")," s\n") #top fin chrono

tmps.ini=Sys.time()         #top chrono
ABC_rejection(.binary_model("./parthy"),prior_matrix,1000,TRUE)
cat("temps mis : ",difftime(Sys.time(),tmps.ini,units="sec")," s\n") #top fin chrono

