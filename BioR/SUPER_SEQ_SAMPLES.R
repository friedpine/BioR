library(RODBC)
GEO_samples <- read.delim("E:/bdyun/00.lab/31.cRNA/samples/GEO_samples_2.txt",header=F)
GEO_samples$names = gsub(".R1.fastq.gz","",GEO_samples$V4)
GEO_fq1 = subset(GEO_samples,V2=="fq1" | V2=="raw_fq1")


library(RODBC)
channel <- odbcConnect("p00", uid="root", pwd="123456")

samples = sqlQuery(channel,"select * from p01.samples")
merged = merge(samples,GEO_fq1,by.x="exp",by.y="V1")
colnames(merged) = c(colnames(merged)[1:12],"GEO_name") 
sqlUpdate(channel,merged[,c(1:8,13)],tablename="p01.samples")

samples = sqlQuery(channel,"select * from p01.samples where length(geo_name)>2")

#HEK_CELLS
FPKM_super_HEK = sqlQuery(channel,"select * from p01.HEK_FPKM_gencode")
for (i in colnames(FPKM_super_HEK)[4:21]){
  geo_name = subset(samples,exp==i)[1,9]
  tmp = FPKM_super_HEK[,c("gene","cat","genename",i)]
  colnames(tmp) = c("ENSEMBL_ID","category","gene_name","FPKM")
  write.table(x=tmp,file=paste(geo_name,"_FPKM_gencode.txt",sep=""))
}

FPKM_super_mouse = sqlQuery(channel,"select * from p01.p01_FPKM")
for (i in colnames(FPKM_super_mouse)[5:56]){
  geo_name = subset(samples,exp==i)[1,9]
  tmp = FPKM_super_mouse[,c("genename",i)]
  colnames(tmp) = c("gene_name","FPKM")
  write.table(x=tmp,file=paste(geo_name,"_FPKM_gencode.txt",sep=""))
}

