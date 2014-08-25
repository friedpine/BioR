library(data.table)
library(RODBC)

SNV_finder = function(channel,samples,control_id,control_type,tablename,target_tables,depth,target_depth){
  controls = samples[control_id]
  targets = samples[-control_id]
  samples_DP = paste(",",samples,"_DP",sep="",collapse="")
  sql_head = paste("select chr,pos,Ref,Alt,Qual,DP,DP_alt",samples_DP," from",tablename,"where",sep=" ")    
  cond_depth = paste(samples,"_DP>=",depth," and ",sep="",collapse = '')
  cond_control_type = paste(controls,"_G = '",control_type,"' and ",sep="",collapse = '')
  cond_control_qual = paste("(",controls,"_GQ>95 or ",controls,"_GQ>3*",controls,"_DP) and ",sep="",collapse = '')
  sql_command = gsub(" and $",";",paste(sql_head,cond_depth,cond_control_type,cond_control_qual))
  SNV_Depth_Ctrl_Type = sqlQuery(channel,sql_command)
  SNV_Depth_Ctrl_Type$POS = paste(as.character(SNV_Depth_Ctrl_Type$chr),SNV_Depth_Ctrl_Type$pos,sep=":")
  INFO_Targets = list()
  
  
  return(SNV_Depth_Ctrl_Type)
}

