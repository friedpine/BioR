library(data.table)
library(RODBC)
library(plyr)
library(gridExtra)
library(reshape2)

seq_base_frequency = function(seqlists){
  apa_seq_mat = t(sapply(strsplit(as.character(seqlists),""),as.character))
  apa_seq_ATCG_counts = sapply(c('A','T','C','G'),function(x) colSums(apa_seq_mat==x)/length(seqlists))  
  return(apa_seq_ATCG_counts)
}

motif_frequency = function(seqlists,motifs){
  seqlists = as.character(seqlists)
  out = sapply(motifs,function(motif) as.numeric(lapply(seqlists,function(x) gregexpr(motif,x)[[1]][[1]])))
  return(out)
}

ATCG_MOTIF_FREQ = function(apa_cluster,filename){
  apa_seq_ATCG_counts = seq_base_frequency(apa_cluster$flankseq200)
  motif_frq = motif_frequency(apa_cluster$flankseq200,motif_rara)
  pdf(filename,width=10,height=15)
  p1 = ggplot(melt(apa_seq_ATCG_counts),aes(x=Var1,y=value,color=Var2))+
    geom_line(size=1)+
    ylim(0,0.6)
  p2 = ggplot(melt(motif_frq),aes(value,color=Var2))+
    geom_density(size=1)+
    xlim(1,200)+ylim(0,0.08)
  print(grid.arrange(p1,p2,ncol=1))
  dev.off()
}

APAs_annotated_cluster_pos = function(channel,aps_cluster_table,utr3_distinct_table,range){
  sql = paste("SELECT a.* from ",aps_cluster_table," a join ",utr3_distinct_table,
    " b on a.gene = b.gene AND a.pos+",range,">b.pos and a.pos-",range,"<b.pos",sep="")
  apa_cluster_annotated_id = sqlQuery(channel,sql)
  return(apa_cluster_annotated_id)
}

APAs_unannotated_motif_pos = function(channel,apas_cluster_table,anno_id,motifs,range){
  apa_cluster_all = sqlQuery(channel,paste("select * from",apas_cluster_table))
  apa_cluster_unannotated = subset(apa_cluster_all,!(id %in% anno_id))
  temp = melt(motif_frequency(apa_cluster_unannotated$flankseq200,motifs))
  apa_cluster_unannotated_motif50_100 = apa_cluster_unannotated[unique(subset(temp,value<range[2] & value>range[1])$Var1),]
  return(apa_cluster_unannotated_motif50_100)
}

Genes_with_multiple_APAs = function(channel,apas,table){
  dt_apas = data.table(apas) 
  group_by_gene = dt_apas[,list(gene=gene[1],chr=chr[1],strand=strand[1],counts=length(pos),apa1=min(pos),apa2=max(pos)),by=gene]
  selected = group_by_gene[counts>1 & apa2-apa1>150]
  out = selected[,c(1,3,4,6,7),with=FALSE]
  sqlSave(channel,dat=out,tablename=table,append=F,rownames=F,colnames=F)
}

Gene_multi_heatmap = function(channel,tablename,near_col,dist_col,exp_cut,sample_low,filename){
  library(gplots)
  data = sqlQuery(channel,paste("select * from",tablename))
  data_sum = data[,dist_col]+data[,near_col]
  data_sum_cut = rowSums(data_sum>=exp_cut)>=sample_low
  distal_by_near = data[,dist_col]/(data[,near_col]+1)
#  distal_by_near = data[,near_col]/(data[,dist_col]+1)
  print(distal_by_near[1:10,])
  distal_by_near_log10 = log10(distal_by_near+0.01)
  heatdata = distal_by_near_log10[data_sum_cut,]
  print(heatdata)
  pdf(filename)
  #dd <- as.dendrogram(hclust(as.dist((1 - cor(t(heatdata)))/2)))
  #heatmap.2(as.matrix(heatdata),dendrogram="row",scale="row",density.info='none',Colv=FALSE,Rowv=dd,col=bluered(100),trace='none')
  heatmap.2(as.matrix(heatdata),dendrogram="row",scale="row",density.info='none',Colv=FALSE,col=bluered(100),trace='none')
  dev.off()
}

Gene_multi_heatmap2 = function(channel,data,near_col,dist_col,exp_cut,sample_low,filename){
  library(gplots)
  data_sum = data[,dist_col]+data[,near_col]
  data_sum_cut = rowSums(data_sum>=exp_cut)>=sample_low
  distal_by_near = data[,dist_col]/(data[,near_col]+1)
  print(distal_by_near[1:10,])
  distal_by_near_log10 = log10(distal_by_near+0.01)
  heatdata = distal_by_near_log10[data_sum_cut,]
  pdf(filename)
  #dd <- as.dendrogram(hclust(as.dist((1 - cor(t(heatdata)))/2)))
  #heatmap.2(as.matrix(heatdata),dendrogram="row",scale="row",density.info='none',Colv=FALSE,Rowv=dd,col=bluered(100),trace='none')
  heatmap.2(as.matrix(heatdata),dendrogram="row",scale="row",density.info='none',Colv=FALSE,col=bluered(100),trace='none')
  dev.off()
}

Gene_multi_density_by_species = function(channel,data,near_col,dist_col,sample_names,sample_table,low_cut,high_cut,filename){
  library(ggplot2)
  library(reshape2)
  samples = sqlQuery(channel,paste("select sample,type from",sample_table))
  data_sum = data[,dist_col]+data[,near_col]
  distal_by_near = data[,dist_col]/(data[,near_col]+1)
  colnames(data_sum) = sample_names
  colnames(distal_by_near) = sample_names
  melt_data_sum = melt(data_sum)
  melt_distal_by_near = melt(distal_by_near)
  melt_dbn = melt_distal_by_near[melt_data_sum$value>low_cut,]
  melt_dbn = merge(melt_dbn,samples,by.x="variable",by.y="sample")
  melt_dbn$log_data = pmax(log10(melt_dbn$value),-2)
  return(melt_dbn)
#  ggplot(melt_dbn,aes())+
#    geom_density()+
#    xlim(-2,4)
}







