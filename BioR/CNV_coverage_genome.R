library(ggplot2)
library(gridExtra)

plot_genomes_single = function(datas,xmax,ymax,filename){
  pdf(file=filename,width=10,height=15)
  print(
    ggplot(datas, aes(bin_id,counts_p,color=chr)) +
      geom_line() +
      geom_abline(intercept=mean(datas$counts_p),slope=0)+
      facet_grid(chr ~ .) +
      ggtitle("NO_AMPLIFICATION_CNV")+
      ylim(0,ymax))
  #USING_POINTS
  print(
    ggplot(datas, aes(bin_id,counts_p,color=chr)) +
      geom_point(size=1) +
      geom_abline(intercept=mean(datas$counts_p),slope=0)+
      facet_grid(chr ~ .) +
      ggtitle("NO_AMPLIFICATION_CNV")+
      ylim(0,ymax))
  #USING_HIST
  print(
    ggplot(datas, aes(counts,fill=chr)) +
      geom_histogram(binwidth = xmax/100) +
      facet_grid(chr ~ .) +
      ggtitle("NO_AMPLIFICATION_CNV")+
      xlim(0,xmax))
  dev.off()
}