library(RODBC)
library(gplots)
library(ggplot2)
library(plyr)
library(matrixStats)


nonzero_std_mean = function(x,low_cut,low_cell_number,rec){
  x = as.matrix(x)
  data_nonzero = x[rowSums(x>=low_cut)>=low_cell_number,]
  df_out = data.frame(Mean=rowMeans(data_nonzero),Sd=rowSds(data_nonzero))
  df_out$log10 = log10(df_out$Mean)
  df_out$CV = df_out$Sd/df_out$Mean
  df_out$rec=rec
  return(df_out)
}

get_3_sigma_line = function(data,xmin,xmax,segs,rec){
  flank_length = (xmax-xmin)/(2*segs)
  poses = (1:segs*2-1)*flank_length+xmin
  df = data.frame(id=1:segs,poses=poses)
  data$group=as.integer((data$log10-poses[1]+flank_length)/(flank_length*2))+1
  CV_M_std = ddply(data,c("group"),summarize,CV_mean=mean(CV),CV_sd=sd(CV))
  CV_M_std = subset(CV_M_std,group %in% df$id)
  df$CV_3std = CV_M_std$CV_mean+3*CV_M_std$CV_sd
  model.res = glm(formula = CV_3std ~ poses, family = poisson, data = df)
  df$esti = model.res$fitted.values
  df$rec = rec
  return(df)
}

random_group_cells_get_3sigma_line = function(data,group_size,group_count){
  total_samples = length(data[1,])
  sample_info = sapply(1:group_count,function(x) sample(1:total_samples,group_size))
  transform_mat = matrix(0,total_samples,group_count)
  for (i in 1:group_count){
    for (j in sample_info[,i]){
      transform_mat[j,i]=1/group_size
    }
  }
  data_group = as.matrix(data) %*% transform_mat
  return(get_3_sigma_line(nonzero_std_mean(data_group,'8pg'),0,3.6,18))
}

plot_scatter_3sigma = function(datas,counts,recs,ymax,filename){
  CVs = lapply(1:counts,function(x) nonzero_std_mean(datas[[x]],1,2,recs[x]))
  CVs_stat = lapply(1:counts,function(x) get_3_sigma_line(CVs[[x]],0,3.6,18,recs[x]))
  df_CVs_stat = ldply(CVs_stat,data.frame)
  df_CVs = ldply(CVs,data.frame)
  pdf(filename,width=15,height=15)
  print(ggplot(df_CVs,aes(log10,CV))+
    geom_point(aes(color=rec))+
    geom_line(data=df_CVs_stat,aes(poses,CV_3std,size=1,color=rec))+
    ylim(0,ymax))
  dev.off()
#  return(list(df_CVs,df_CVs_stat))
}

