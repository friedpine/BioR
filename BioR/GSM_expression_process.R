combine_columes = function(samples,prefix,colume_id,type){
  cmd = sapply(samples$sample,function(x) paste("less ",x,"_* | awk '{print $",colume_id,"}' >",prefix,x,".txt",sep=""))
  cmd = c(cmd,paste("paste ",prefix,"* >",prefix,type,".txt",sep=""))
  write(cmd,"cmd.txt")
}

