library(matrixStats)
library(ggplot2)
library(data.table)
library(RODBC)

#USING FU DATA
hist(test$V2)
df_2K = subset(test,V2<3000&V2>2000)

df_plot = data.frame(Mean=colMeans(df_2K[3:22]),
                     Std=colSds(as.matrix(df_2K[3:22])))
write.table(df_2K,"df_2K.txt")

#SAMTOOLS DEPTH FILE
dt_depth = data.table(depth)
colnames(dt_depth) = c("transc","pos","depth")

channel <- odbcConnect("my", uid="root", pwd="123456")
refseq = sqlQuery(channel,"select * from mm10_enseall.Refseq")
refseq = refseq[!duplicated(refseq$transc),]

dt_transc = dt_depth[,list(cc=sum(depth)),by=transc]
dt_transc = subset(dt_transc,transc %in% refseq$transc)
dt_transc_info = merge(dt_transc,refseq,by="transc")
dt_transc_info = subset(dt_transc_info,cc>=length*10 & isoform==1)

dt_depth_select = dt_depth[transc %in% dt_transc_info$transc]
dt_depth_select = merge(dt_depth_select,dt_transc_info,by="transc")
dt_depth_select$binID = as.integer(dt_depth_select$pos*20/dt_depth_select$length)+1

dt_bined = dt_depth_select[,list(bin_DP=sum(depth)/mean(cc),
                                 length=mean(length)),by=list(transc,binID)]
dt_bined_2K = subset(dt_bined,length<3000 & length>2000)
dt_bined_2K_plot = dt_bined_2K[,list(M=mean(bin_DP),
                                     S=sd(c(rep(0,1286-length(bin_DP)),bin_DP))),
                               by=list(binID)]
ggplot(dt_bined_2K_plot,aes(x=binID,y=M))+
  geom_line()+
  geom_errorbar(aes(ymin=M-S,ymax=M+S))


divide20 = function(df_data,length,total){
  df_data$V4=as.integer(df_data$V2*20/(length+20)+1)
  return(df_data)
}

divide20(dt_depth[V1=="NM_001001806"],3535,293744)
