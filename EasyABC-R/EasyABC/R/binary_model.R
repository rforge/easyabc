## model wrapper
################
binary_model<-function(command) {
  invoke<-function(param) {
    write.table(param,file="input",row.names=F,col.names=F,quote=F)
    system(command)
    read.table("output",h=F)
  }
}
