library(data.table)
library(RODBC)
library(plyr)

SNV_finder = function(channel,samples,control_id,control_type,tablename,target_tables,depth,target_depth,alt_ratio){
  controls = samples[control_id]
  targets = samples[-control_id]
  samples_DP = paste(",",samples,"_DP",sep="",collapse="")
  sql_head = paste("select chr,pos,Ref,Alt,Qual,DP,DP_alt",samples_DP," from",tablename,"where",sep=" ")    
  cond_depth = paste(samples,"_DP>=",depth," and ",sep="",collapse = '')
  cond_control_type = paste(controls,"_G = '",control_type,"' and ",sep="",collapse = '')
  cond_control_qual = paste("(",controls,"_GQ>95 or ",controls,"_GQ>3*",controls,"_DP) and ",sep="",collapse = '')
  sql_command = gsub(" and $",";",paste(sql_head,cond_depth,cond_control_type,cond_control_qual))
  print(sql_command)
  SNV_Depth_Ctrl_Type = sqlQuery(channel,sql_command)
  SNV_Depth_Ctrl_Type$POS = paste(as.character(SNV_Depth_Ctrl_Type$chr),SNV_Depth_Ctrl_Type$pos,sep=":")
  INFO_Targets = list()
  sql_target = function(char1,int2,int3){paste("select chr,pos,DP,DP_alt from ",char1," where DP>=",int2," and DP_alt>=",int3," and DP_alt>=DP*",alt_ratio,sep="")}
  INFO_Targets = lapply(target_tables,function(x){sqlQuery(channel,sql_target(x,depth,target_depth))})
  INFO_Targets = llply(INFO_Targets,function(x) cbind(DP=x$DP,DP_alt=x$DP_alt,POS=paste(x$chr,x$pos,sep=":")))
  INFO_Targets[[length(INFO_Targets)+1]] = SNV_Depth_Ctrl_Type
  merged = Reduce(function(...) merge(...,by="POS",all=F),INFO_Targets)
  merged_info = llply(INFO_Targets,function(x) length(merge(x,INFO_Targets[[length(INFO_Targets)]],by="POS")$POS))
  names(merged_info) = c(targets,"SNV_Depth_Ctrl_Type")
  infos = append(list(depth=depth,ratio=alt_ratio,target_depth=target_depth,overlap=length(merged$POS)),merged_info)
  return(list(infos,merged))
}


SNV_Double_Hit = function(channel,samples,tablename,target_tables,depth,target_depth,alt_ratios,ratio_rec){
  samples_DP = paste(",",samples,"_DP",sep="",collapse="")
  sql_head = paste("select chr,pos,Ref,Alt,Qual,DP,DP_alt",samples_DP," from",tablename,"where",sep=" ")    
  cond_depth = paste(samples,"_DP>=",depth," and ",sep="",collapse = '')
  sql_command = gsub(" and $",";",paste(sql_head,cond_depth))
  print(sql_command)
  SNV_Depth_Ctrl_Type = sqlQuery(channel,sql_command)
  SNV_Depth_Ctrl_Type$POS = paste(as.character(SNV_Depth_Ctrl_Type$chr),SNV_Depth_Ctrl_Type$pos,sep=":")
  INFO_Targets = list()
  sql_target = function(char1,sample,int2,int3,range){a = paste("select chr,pos,DP as ",sample,"_DP,DP_alt as ",sample,"_DP_alt from ",char1," where DP>=",int2," and DP_alt>=",int3," and DP_alt>=DP*",range[1]," and DP_alt<=DP*",range[2],sep="");print(a)}
  INFO_Targets = lapply(1:length(samples),function(x){sqlQuery(channel,sql_target(target_tables[x],samples[x],depth,target_depth,alt_ratios[[x]]))})
  INFO_Targets = llply(INFO_Targets,function(x) cbind(DP=x$DP,DP_alt=x$DP_alt,POS=paste(x$chr,x$pos,sep=":")))
  INFO_Targets[[length(INFO_Targets)+1]] = SNV_Depth_Ctrl_Type
  merged = Reduce(function(...) merge(...,by="POS",all=F),INFO_Targets)
  merged_info = llply(INFO_Targets,function(x) length(merge(x,INFO_Targets[[length(INFO_Targets)]],by="POS")$POS))
  names(merged_info) = c(samples,"SNV_Depth_Ctrl_Type")
  infos = append(list(depth=depth,ratio=ratio_rec,target_depth=target_depth,overlap=length(merged$POS)),merged_info)
  return(list(infos,merged))
}

SNV_Tested_by_Other_Method = function(SNVs,tables_names){
  
}




