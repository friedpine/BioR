library(ggplot2)
library(reshape2)
library(data.table)
library(DNAcopy)
source('~/Documents/BioR/BioR/CNV_coverage_genome.R')

agg = function(x){sapply(1:length(x),function(i){sum(x[1:i])})}

norm_by_samples = function(depth,samples,refs){
  samples_all = c(samples,refs)
  samples_all = samples_all[!duplicated(samples_all)]
  df_norm = as.data.frame(apply(depth[,samples_all],2,function(x){x*1000000/sum(x)}))
  print(colnames(df_norm))
  print(refs)
  df_norm$ref = rowMeans(df_norm[,refs])
  reserved_pos = df_norm$ref<=20
  df_cnv = as.data.frame(matrix(2,nrow = dim(df_norm)[1],ncol=length(samples)))
  colnames(df_cnv) = samples
  df_cnv[!reserved_pos,] = as.data.frame(apply(df_norm[,samples],2,function(x){2*x/df_norm$ref}))[!reserved_pos,]
  df_cnv$x = 1:dim(df_norm)[1]
  return(df_cnv)
}

alter_CNV_size = function(depth,samples,size){
  depth = data.table(depth);
  depth$bin_id = as.integer(depth$bin_id/size);
  depth=depth[, lapply(.SD, mean), by = list(chr,bin_id)];
  depth = as.data.frame(depth);
  depth$pos = 1:length(depth$bin_id)
  return(depth)
}

CNVs_segs = function(CNV_all,samples,sample_norm,size,sexfactor){
  df_CNVs = norm_by_samples(CNV_all,samples,sample_norm)
  df_CNVs = alter_CNV_size(cbind(CNV_all[,1:2],df_CNVs),samples,size)
  
  df_CNVs[df_CNVs$chr=='chrX',samples] = df_CNVs[df_CNVs$chr=='chrX',samples]*(sexfactor[1]/2);
  df_CNVs[df_CNVs$chr=='chrY',samples] = df_CNVs[df_CNVs$chr=='chrY',samples]*(sexfactor[2]/2);
  
  df_CNVs[,samples] = as.matrix(df_CNVs[,samples])%*%diag(2/apply(df_CNVs[,samples],2,median))
  
  CNA.object = CNA(log(df_CNVs[,samples]+0.01),
                   df_CNVs$chr,df_CNVs$bin_id,
                   data.type="logratio",sampleid=samples)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  smoothed.CNA.object[,samples] = exp(smoothed.CNA.object[,samples])-0.01
  ssc <- segment(smoothed.CNA.object, verbose=1)
  ssc$output$seg.mean[ssc$output$seg.mean>=4]=4
  
  agg = function(x){sapply(1:length(x),function(i){sum(x[1:i])})}
  
  segs_plot = ddply(ssc$output,.(ID),transform,end=agg(num.mark))
  colnames(segs_plot)[1] = 'variable'
  
  plot_CNVs_segs(df_CNVs,segs_plot,samples,1)
  return(list(cnvs=df_CNVs,seg=segs_plot))
}

plot_CNVs = function(depth,samples,size){
  depth = data.table(depth);
  depth$bin_id = as.integer(depth$bin_id/size);
  depth=depth[, lapply(.SD, mean), by = list(chr,bin_id)];
  depth = as.data.frame(depth);
  depth$pos = 1:length(depth$bin_id);
  
  df_plot = melt(depth[,c("pos",samples)],"pos")
  
  p1 = ggplot(df_plot,aes(x=pos,y=value))+
    geom_point(size=1)+
    facet_grid(variable ~ .) +
    geom_abline(intercept = 2, slope = 0)+
    xlim(1,max(df_plot$pos))+
    ylim(0,4)+
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_rect(fill=NA,color="black", size=0.5,
                                  linetype="solid"))+
    scale_x_continuous(expand = c(0,0))
  print(p1)
  return(depth)
}

plot_CNVs_segs = function(depth,segs,samples,size){
  
  depth = data.table(depth);
  print(depth[1]);
  depth$bin_id = as.integer(depth$bin_id/size);
  depth=depth[, lapply(.SD, mean), by = list(chr,bin_id)];
  depth = as.data.frame(depth);
  
  depth$pos = 1:length(depth$bin_id)
  
  df_plot = melt(depth[,c("pos",samples)],"pos")
  
  df_plot$value[df_plot$value>=4] = 3.99
  
  df_rect = ddply(depth[,c("chr","pos")],.(chr), summarise, mm=min(pos), ma=max(pos));
  df_rect$c = as.factor((1:dim(df_rect)[1] %%2));
  print(df_rect)
  print(segs)
  
  #labels
  df_label = ddply(depth[,c("chr","bin_id","pos")],.(chr),summarise,label=strsplit(as.character(chr[1]),"chr")[[1]][2],pos=as.integer(median(pos)))
  df_label = df_label[seq(1,length(df_rect[,1]),2),];
  
  p1 = ggplot(df_plot,aes(x=pos,y=value));
  
  for (x in seq(1,length(df_rect[,1]),2)){
    p1 = p1+geom_rect(xmin=df_rect[x,'mm'],xmax=df_rect[x,'ma'],ymin=0,ymax=4,fill='gray83',alpha=0.1)
  }
  
  p1 = p1+
    geom_point(size=1)+
    geom_abline(intercept = 2, slope = 0)+
    geom_segment(data=segs,aes(x=end-num.mark+1,xend=end-1,y=seg.mean,yend=seg.mean),size=1,colour="red")+
    facet_grid(variable ~ .) +
    scale_x_continuous(limits=c(1,max(df_plot$pos)),expand = c(0,0),breaks=df_label$pos,labels=df_label$label)+
    ylim(0,4)+
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text.x = element_text(size=16, angle=0, color="black"),
      panel.border = element_rect(fill=NA,color="black", size=0.5,
                                  linetype="solid"))+
    xlab("")+
    ylab("Copy Number")
  
    print(p1)
}

