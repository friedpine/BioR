library(sqldf)

CNV_Segs_and_HeterRatio = function(segs,binlen,sample,ratioPath,minSites){
  df_r = subset(segs,variable==sample)[,c(2,3,4,6)]
  colnames(df_r) = c("chr","start","end","CNV")
  
  hp <- read.delim(ratioPath, header=FALSE)
  hp$chr = sapply(strsplit(as.character(hp$V3),"_"),"[",1)
  hp$pos = sapply(strsplit(as.character(hp$V3),"_"),"[",2)
  hp$binid = as.integer(as.integer(hp$pos)/binlen)
  hp$ratio = as.double(sapply(strsplit(as.character(hp$V5),"#"),"[",2))
  hp$ratioABS = 0.5-abs(hp$ratio-0.5)
  hp = hp[,6:10]
  joined = sqldf("select * from df_r a join hp b on a.chr=b.chr and b.binid>=a.start and b.binid<=a.end")

  tt = ddply(joined,.(chr,start,end),summarise,CNV=CNV[1],med=median(ratio),median_ratio = median(ratioABS),counts=length(ratioABS))
  tt$sample = sample
  return(tt)
}

chrs = c("chr1", "chr2",  "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
         "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20", "chr21", "chr22", "chrX", "chrY")

ratioPath = "~/Documents/project/LH/DATA_CNVs/HerterRatio/Haplo_YJY_T68.txt"
binlen = 3000000

#YJY
path = "~/Documents/project/LH/DATA_CNVs/HerterRatio/Haplo_YJY_"
samples = c("MT68","MT68_2","MT69","MT69_2","MT70","MT70_2","MT71",
            "MT71_2","MT72","MT72_2","MT73","MT73_2")
samples2 = c("T68","T68_2","T69","T69_2","T70","T70_2","T71",
            "T71_2","T72","T72_2","T73","T73_2")
df1 = out8$seg

#ZWJ
path = "~/Documents/project/LH/DATA_CNVs/HerterRatio/Haplo_ZWJ_"
samples=c("MT16","MT17","MT3_M")
samples2=c("T16","T17","T3")
df1 = out7$seg

#ZCG
path = "~/Documents/project/LH/DATA_CNVs/HerterRatio/Haplo_ZCG_"
samples=c("MT20","MT21","MT22","MT20a","MT21a","MT20a_2","MT21a_2")
samples2=c("T20","T21","T22","T20a","T21a","T20a_2","T21a_2")
df1 = out5$seg

#ZES
path = "~/Documents/project/LH/DATA_CNVs/HerterRatio/Haplo_ZES_"
samples=c("MT6","MT7")
samples2=c("T6","T7")
df1 = out6$seg

#LXY
path = "~/Documents/project/LH/DATA_CNVs/HerterRatio/Haplo_LXY_"
samples=c("MT30","MT31","MT32","MT33")
samples2=c("T30","T31","T32","T33")
samples=c("MT30","MT31","MT32")
samples2=c("T30","T31","T32")
df1 = out1$seg

#XLY
path = "~/Documents/project/LH/DATA_CNVs/HerterRatio/Haplo_XLY_"
samples=c("MT18","MT19")
samples2=c("T18","T19")
df1 = out3$seg

#YZY
path = "~/Documents/project/LH/DATA_CNVs/HerterRatio/Haplo_YZY_"
samples=c("MT4","MT5")
samples2=c("T4","T5")
df1 = out4$seg

####
res = list()
for(x in 1:length(samples)){
  sample = samples[x]
  res[[x]] = CNV_Segs_and_HeterRatio(df1,3000000,sample,
                                     paste(path,samples2[x],".txt",sep=""),50)
}
dfm = Reduce(function(...) merge(..., all=T), res)
dfm$ten = dfm$median_ratio*10

ggplot(subset(dfm,counts>=10),aes(CNV,median_ratio,color=sample))+
  geom_point()+
  scale_x_continuous(limits=c(0,4))

#Kmeans!!!
dfm$cluster = as.factor(kmeans(dfm[,c("CNV","ten")],3,iter.max=100)$cluster)
ggplot(subset(dfm,counts>=10),aes(CNV,median_ratio,color=cluster))+
  geom_point()+
  scale_x_continuous(limits=c(0,4))

#Plot!!!
df_plot = merge(dfm,df_pos_3M,by.x="chr",by.y="chrom")[,c("sample","chr","agg","start","end","CNV","cluster")]
df_plot$start2 = df_plot$start+df_plot$agg
df_plot$end2 = df_plot$end+df_plot$agg
df_plot$cluster = as.factor(df_plot$cluster)
ggplot(df_plot,aes(xmin=start2,xmax=end2,ymin=0,ymax=4,fill=cluster))+
  facet_grid(sample ~ .)+
  geom_rect()+
  scale_x_continuous(breaks=c(0,200,400,600,800,1000))









  



