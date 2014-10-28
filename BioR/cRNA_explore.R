library(data.table)
library(ggplot2)
library(reshape2)
library(gplots)
library(plyr)

cRNA_counts_and_gene_FPKM = function(channel,event_table,FPKM_table,sample_table,g2t_table){
	sql1 = paste("select a.*,b.exp,c.gene from ",event_table," a join ",sample_table," b join ",g2t_table," c on a.sample=b.sample and a.transc=c.transc",sep="")
	sql2 = paste("select * from ",FPKM_table)
	df_event = sqlQuery(channel,sql1)
	df_FPKM = sqlQuery(channel,sql2)
	print(df_event[1:5,])
	print(length(df_event[,1]))
	df_FPKM = subset(df_FPKM,gene %in% df_event$gene)
	df_FPKM = df_FPKM[,!(names(df_FPKM) %in%c("genename","cate","LENGTH"))]
	df_FPKM_melt = melt(df_FPKM)
	colnames(df_FPKM_melt) = c("gene","exp","FPKM")
	print(df_FPKM_melt[1:5,])
	df_event_FPKM = merge(df_event,df_FPKM_melt,by=c("gene","exp"),all.x=T)
	#return(df_event_FPKM)
	dt_event_FPKM = data.table(df_event_FPKM)
	dt_event_summ = dt_event_FPKM[,list(crna_reads=length(event),FPKM=FPKM[1]),by=c("gene","sample")]
	return(dt_event_summ)
}

cRNA_sample_pairs_overlap = function(channel,event_table,samples1,samples2,filename){
	sql1 = paste("select sample,event from ",event_table," where t_dis=0")
	df_event = sqlQuery(channel,sql1)
	df_event$c =1
	print(df_event[1:10,])
	pdf(filename)
	for(x in 1:length(samples1)){
		sample1 = samples1[x]
		sample2 = samples2[x]
		print(sample1)
		circs1 = unique(subset(df_event,sample==sample1)[,2])
		circs2 = unique(subset(df_event,sample==sample2)[,2])
		venn(list(S1=circs1,S2=circs2))
	}
	dev.off()
	return(ddply(df_event,.(sample,event),summarize,sample=sample[1],event=event[1],count=sum(c)))
}



density_of_repeats = function(datalist,datanames,repeat_type){
	df = data.frame(count=c(),len=c(),name=c())
	for(x in 1:length(datanames)){
		dt = datalist[[x]]
		sample = datanames[x]
		result = data.frame(dt[,list(count=sum(family==repeat_type),
                   len=intron_end[1]-intron_start[1]),
             by=c("transc","in_id")])
		print(result[1:10,])
		print(sample)
		result$typename = sample
		df = rbind(df,result[,3:5])} 
	df$density = df$count*1000/df$len
	return(df)
}

repeat_pairing_of_circRNA = function(datalist,repeat_type,total_circ){
	dt_pairings = data[,list(),by=c("event")]		
}

