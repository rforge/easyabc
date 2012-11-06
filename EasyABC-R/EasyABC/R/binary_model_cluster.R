## model wrapper
################
binary_model_cluster<-function(command) {
  invoke<-function(param) {
    n_clust=param[1]
    numclust=1+((param[2]-1)%%n_clust)
    nparam=length(param)
    input_file_name=paste("input",numclust,sep="")
    output_file_name=paste("output",numclust,sep="")
    write.table(param[2:nparam],file=input_file_name,row.names=F,col.names=F,quote=F)
    system2(command,args=as.character(numclust))
    file.remove(input_file_name)
    read.table(output_file_name,h=F)
    file.remove(output_file_name)
  }
}
