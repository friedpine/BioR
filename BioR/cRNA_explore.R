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

intron_pair_blast_process = function(channel,introns_table,blast_table,out_table){
	infos = sqlQuery(channel,paste("select * from ",introns_table,sep=""))
	blast = sqlQuery(channel,paste("select * from ",blast_table,sep=""))
	blast_infos = merge(blast,infos,by="event")
	blast_infos$rc_in1_s = blast_infos$s_end+blast_infos$in1_s 
	blast_infos$rc_in1_e = blast_infos$s_start+blast_infos$in1_s
	blast_infos$rc_in2_s = blast_infos$in2_s+blast_infos$q_start
	blast_infos$rc_in2_e = blast_infos$in2_s+blast_infos$q_end
	sqlSave(channel,blast_infos[,c(1,4,5,6,7,12,13,14,15,18,23:26)],tablename=out_table,rownames=F,colnames=F)
}

#(channel,'mm_crna_explore.b3_rc_poses','mm_crna_explore.02_circ_left_intron_info_repeatmask','mm_crna_explore.02_circ_right_intron_info_repeatmask','mm_crna_explore.b4_rc_repeatmasks')
intron_pair_blsat_to_repeat = function(channel,rc_table,in1_repeat,in2_repeat,out_table){
	infos1 = sqlQuery(channel,paste("SELECT a.id,a.event,b.family as up_family,b.strand as up_strand FROM ",rc_table," a LEFT JOIN ",in1_repeat," b ON a.transc=b.transc AND a.rc_in1_s>b.start-25 AND a.rc_in1_e<b.end+25 GROUP BY id"))
	infos2 = sqlQuery(channel,paste("SELECT a.id,b.family as down_family,b.strand as down_strand FROM ",rc_table," a LEFT JOIN ",in2_repeat," b ON a.transc=b.transc AND a.rc_in2_s>b.start-25 AND a.rc_in2_e<b.end+25 GROUP BY id"))
	info_m = merge(infos1,infos2,by="id")
	sqlSave(channel,info_m,tablename=out_table,rownames=F,colnames=F)
	}

intron_pair_ReverseComp_summ = function(channel,event_table,rc_blast_table,out_table){
	sqlcmd1 = sprintf("select a.event,b.* from %s a left join %s b on a.event=b.event",event_table,rc_blast_table)
	infos = sqlQuery(channel,sqlcmd1)
	infos$cates = 'Others'
	infos$cates[is.na(infos$id)] = "NO_RC"
	infos$cates[!is.na(infos$up_family) & !is.na(infos$down_family)] = "RC_repeat"
	infos$cates[is.na(infos$up_family) & is.na(infos$down_family) & !is.na(infos$id)] = "RC_Nonrepeat"
	dt_info  = data.table(infos)
	dt_summ = dt_info[,list(no_rc=sum(cates=="NO_RC"),rc_rep=sum(cates=="RC_repeat"),rc_norep=sum(cates=="RC_Nonrepeat")),by=c("event")]
	sqlSave(channel,dt_summ,tablename=out_table,rownames=F,colnames=F)
}


