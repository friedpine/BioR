library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)

plot_genomes_single = function(datas,xmax,ymin,ymax,filename){
  pdf(file=filename,width=10,height=15)
  print(
    ggplot(datas, aes(bin_id,counts_p,color=chr)) +
      geom_line() +
      geom_abline(intercept=mean(datas$counts_p),slope=0)+
      facet_grid(chr ~ .) +
      ggtitle("NO_AMPLIFICATION_CNV")+
      ylim(ymin,ymax))
  #USING_POINTS
  print(
    ggplot(datas, aes(bin_id,counts_p,color=chr)) +
      geom_point(size=1) +
      geom_abline(intercept=mean(datas$counts_p),slope=0)+
      facet_grid(chr ~ .) +
      ggtitle("NO_AMPLIFICATION_CNV")+
      ylim(ymin,ymax))
  #USING_HIST
  print(
    ggplot(datas, aes(counts,fill=chr)) +
      geom_histogram(binwidth = xmax/100) +
      facet_grid(chr ~ .) +
      ggtitle("NO_AMPLIFICATION_CNV")+
      xlim(0,xmax))
  dev.off()
}

noise_distraction_plot = function(datas,level,filename){
  mean1 = mean(datas[,3]) 
  mean2 = mean(datas[,4])
  datas[,3] = datas[,3]*(1000/mean1)
  datas[,4] = datas[,4]*(1000/mean2)
  datas2 = datas[,1:2]
  datas2$counts = datas[,3]-datas[,4]
  datas2$counts_p = pmin(datas2$counts,level*2)
  plot_genomes_single(datas2,2500,-2000,2000,filename)
}


sample_divided_by_ref_plot = function(datas,ymax,filename){
  mean1 = mean(datas[,3]) 
  mean2 = mean(datas[,4])
  datas[,3] = datas[,3]*(1000/mean1)
  datas[,4] = datas[,4]*(1000/mean2)
  datas2 = datas[,1:2]
  datas2$counts = (datas[,3]+100)/(datas[,4]+100)*2
  datas2$counts_p = pmin(datas2$counts,ymax)
  plot_genomes_single(datas2,ymax,0,ymax,filename)
}

multi_samples_plot = function(channel,samples,table_name,filename_pre){
  for (sample in samples){
    datas = sqlQuery(channel,paste("select chr,bin_id,",sample," as counts from ",table_name,sep=""))
    mean_count = mean(datas$counts)
    datas$counts_p = pmin(datas$counts,mean_count*2)
    plot_genomes_single(datas,mean_count*2,0,mean_count*2,paste(filename_pre,"_",sample,".pdf",sep=""))
  }
}

multi_samples_plot_divided_by_ref = function(channel,samples,ref_samples,table_name,filename_pre){
  for (sample in samples){
    ref_sample = ref_samples[which(samples==sample)]
    datas = sqlQuery(channel,paste("select chr,bin_id,",sample,",",ref_sample," from ",table_name,sep=""))
    sample_divided_by_ref_plot(datas,4,paste(filename_pre,"_",sample,".pdf",sep=""))
  }  
}

normalized_by_depth = function(data,infoids,colids){
  norm_factor = diag(length(data[,1])/colSums(data[,colids])) 
  norm_data = as.matrix(data[,colids]) %*% as.matrix(norm_factor)
  colnames(norm_data) = colnames(data)[colids]
  if(infoids==0){
    return(as.data.frame(norm_data))
  }
  return(cbind(data[,infoids],norm_data))
}

distribution_samples = function(data,colids,xrange){
  dev.new()
  dt_plot = melt(data[,colids])
  ggplot(dt_plot,aes(value,fill=variable,color=variable),size=30)+
  geom_density(alpha=0.2)+
  xlim(xrange[1],xrange[2])
}

normalized_by_ColAverage = function(data,infoids,colids,StandradIds){
  colaverage = rowMeans(data[,StandradIds])
  out = data[,infoids]
  for(x in colids){
    out=cbind(out,(data[,x]+0.01)/(colaverage+0.01))
  }
  colnames(out) = colnames(data)[c(infoids,colids)]
  return(out)
}


average_by_multiSegs = function(data,infoids,colids,counts){
  total = length(data[,1])
  half_count = as.integer(counts/2)
  ranges = (half_count+1):(total-half_count)
  out_data = matrix(,length(ranges),length(colids))
  for (x in ranges){
    out_data[x-half_count,]=colMeans(data[(x-half_count):(x+half_count-1),colids])
  }
  out = cbind(data[ranges,infoids],out_data)
  colnames(out) = colnames(data)[c(infoids,colids)]
  return(out)
}

plot_norm_malbac_data = function(data_in,colum_plot,filename_prefix){
  dev.new()
  filename=paste(filename_prefix,"_",colnames(data_in)[colum_plot],".pdf",sep="")
  pdf(file=filename,width=10,height=15)
  datas = data_in[,c(1,2,colum_plot)]
  colnames(datas) = c("chr","bin_id","depth")
  print(ggplot(datas,aes(bin_id,depth,color=chr)) +
  geom_point(size=1) +
  geom_abline(intercept=1,slope=0)+
  facet_grid(chr ~ .) +
  ylim(0,2))
  dev.off()
}

plot_norm_malbac_datas = function(data_in,colums_plot,filename){
  #dev.new()
  filename=paste(filename,".pdf",sep="")
  pdf(file=filename,width=10,height=15)
  for(x in colums_plot){
  datas = data_in[,c(1,2,x)]
  colnames(datas) = c("chr","bin_id","depth")
  print(ggplot(datas,aes(bin_id,depth,color=chr)) +
  geom_point(size=1) +
  geom_abline(intercept=1,slope=0)+
  facet_grid(chr ~ .) +
  ylim(0,2)+
  labs(title = colnames(data_in)[x]))}
  dev.off()
}

depth_corr_coef_matrix = function(data,low_cut,xlimit,plotfile){
  samples = colnames(data)
  counts = length(samples)
  cormat = matrix(0,nrow=counts,ncol=counts)
  colnames(cormat) = samples
  rownames(cormat) = samples
  pdf(plotfile)
  for (x in seq(1,counts-1)){
    for (y in seq(x+1,counts)){
      data_exp = data[,c(x,y)]
      data_nonzero = data_exp[rowSums(data_exp>low_cut)>0,]
      cormat[x,y] = round(cor(data_nonzero[,1],data_nonzero[,2]),2)
      cormat[y,x] = cormat[x,y]
      title_name = paste(samples[x],samples[y],sep="_")
      print(ggplot(data.frame(x=data_nonzero[,1],y=data_nonzero[,2]),aes(x,y))+
        geom_point(size=1)+
        xlim(xlimit[1],xlimit[2])+ylim(xlimit[1],xlimit[2])+
        labs(title = title_name)
        )
    }
    cormat[x,x] = 1
  }
  cormat[counts,counts] = 1
  dev.off()
  return(cormat)
}

whole_genome_multisamples = function(depth,species,binsize,ymax,controls,bings){
  library(data.table)
  depth = data.table(depth);
  depth$bin_id = as.integer(depth$bin_id/binsize);
  depth=depth[, lapply(.SD, sum), by = list(chr,bin_id)];
  depth = as.data.frame(depth);
 
  out = depth[,c('chr','bin_id',bings)];
  print(out[1,]);
  print(length(out[,1]));
 
  if (controls != 0){
  for(x in 1:length(bings)){
    out[,bings[x]] = 0;
    bing = depth[,bings[x]];
    control = depth[,controls[x]];
    out[,bings[x]] = 2*((bing+10*sum(bing)/sum(control))/(control+10))/(sum(bing)/sum(control));
  }}
 
  for(x in 1:length(bings)){
     temp = out[,bings[x]]
     out[temp>ymax,bings[x]] = ymax
}

  print(out[1,]);
  print(length(out[,1]));
 
  if(species=='hg19'){
    chrs = c("chr1", "chr2",  "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20", "chr21", "chr22", "chrX", "chrY");
  }
 
  if(species=='mm10'){
    chrs = c("chr1", "chr2",  "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chrX", "chrY");
  }
 
  out$chr=factor(out$chr,levels=chrs);
 
  out$x=1:length(out[,1]);
  out_col_length = length(colnames(out));
  print(out_col_length);
 
  #out[out$chr=='chrX',3:(out_col_length-1)]=out[out$chr=='chrX',3:(out_col_length-1)]/2;
  #out[out$chr=='chrY',3:(out_col_length-1)]=out[out$chr=='chrY',3:(out_col_length-1)]/2;
 
  out1 = out[,c(out_col_length,3:(out_col_length-1))]
  df_out = melt(out1,'x');
 
 
  print(out1[1,]);
  print(df_out[1,]);
  print(dim(df_out));
 
  df_rect = ddply(out[,c(1:2,out_col_length)],.(chr),summarise,mm=min(x),ma=max(x));
  df_rect$c = as.factor((1:length(chrs)) %%2);
  print(df_rect);
 
  p1 = ggplot(df_out,aes(x=x,y=value));
 
  for (x in seq(1,length(df_rect[,1]),2)){
    p1 = p1+geom_rect(xmin=df_rect[x,'mm'],xmax=df_rect[x,'ma'],ymin=0,ymax=4,fill='gray83',alpha=0.1)
  }
 
  p1 = p1+geom_point(size=1)+
    facet_grid(variable ~ .) +
    geom_abline(intercept = 2, slope = 0)+
    xlim(1,max(df_out$x))+
    ylim(0,ymax)+
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_rect(fill=NA,color="black", size=0.5,
                                  linetype="solid"));
 
  print(p1);
}





