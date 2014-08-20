library(matrixStats)
library(ggplot2)
library(data.table)
library(RODBC)

#IMPORTING SAMTOOLS DEPTH FILE (samtools depth BWA_MAP2Refseq.bam > depth.txt)
depth <- read.table("E:/PATH/TO/DepthFILE", quote="\"")
dt_depth = data.table(depth)
colnames(dt_depth) = c("transc","pos","depth")

#LOAD REFSEQ INFORMATION: GENENAME,length
channel <- odbcConnect("my", uid="root", pwd="123456")
refseq = sqlQuery(channel,"select * from mm10_enseall.Refseq")
refseq = refseq[!duplicated(refseq$transc),]

#SELECT transcript with mean depth>= 10 and Only one isoform
dt_transc = dt_depth[,list(cc=sum(depth)),by=transc]
dt_transc = subset(dt_transc,transc %in% refseq$transc)
dt_transc_info = merge(dt_transc,refseq,by="transc")
dt_transc_info = subset(dt_transc_info,cc>=length*10 & isoform==1)

#Calculate the relative sequencing depth of each bin from each gene 
dt_depth_select = dt_depth[transc %in% dt_transc_info$transc]
dt_depth_select = merge(dt_depth_select,dt_transc_info,by="transc")
dt_depth_select$binID = as.integer(dt_depth_select$pos*20/dt_depth_select$length)+1
dt_bined = dt_depth_select[,list(bin_DP=sum(depth)/mean(cc),
                                 length=mean(length)),by=list(transc,binID)]

#FUNCTION 函数化
depth2bined = function(dt_depth,refseq,depth_cutoff){
    #SELECT transcript with mean depth>= 10 and Only one isoform
    dt_transc = dt_depth[,list(cc=sum(depth)),by=transc]    
    dt_transc = subset(dt_transc,transc %in% refseq$transc)
    dt_transc_info = merge(dt_transc,refseq,by="transc")
    dt_transc_info = subset(dt_transc_info,cc>=length*depth_cutoff & isoform==1)
    #Calculate the relative sequencing depth of each bin from each gene 
    dt_depth_select = dt_depth[transc %in% dt_transc_info$transc]
    dt_depth_select = merge(dt_depth_select,dt_transc_info,by="transc")
    dt_depth_select$binID = as.integer(dt_depth_select$pos*20/dt_depth_select$length)+1
    dt_bined = dt_depth_select[,list(bin_DP=sum(depth)/mean(cc),
                                 length=mean(length)),by=list(transc,binID
    return(dt_bined)
}


# Prepare ploting data of genes around 2-3K
dt_bined = depth2bined(dt_depth,refseq,10)


dt_bined_2K = subset(dt_bined,length<3000 & length>2000)
dt_bined_2K_plot = dt_bined_2K[,list(M=mean(bin_DP),
                                     S=sd(c(rep(0,1286-length(bin_DP)),bin_DP))),
                               by=list(binID)]
                               
ggplot(dt_bined_2K_plot,aes(x=binID,y=M))+
  geom_line()+
  geom_errorbar(aes(ymin=M-S,ymax=M+S))
