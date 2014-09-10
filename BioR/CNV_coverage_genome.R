library(ggplot2)
library(gridExtra)

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




