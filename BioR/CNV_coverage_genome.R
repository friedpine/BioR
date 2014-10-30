library(ggplot2)
library(gridExtra)
library(reshape2)

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

plot_data = function(data_in,colum_plot,filename_prefix){
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




