## model wrapper
################
binary_model_cluster<-function(command) {
  invoke<-function(param) {
    n_clust=param[1]
    numclust=1+((param[2]-1)%%n_clust)
    nparam=length(param)
    write.table(param[2:nparam],file=paste("input",numclust,sep=""),row.names=F,col.names=F,quote=F)
    system2(command,args=as.character(numclust))
    read.table(paste("output",numclust,sep=""),h=F)
  }
}
