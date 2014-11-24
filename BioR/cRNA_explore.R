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
	print(sql1)
	print(sql2)
	print(df_event[,])
	print(length(df_event[,1]))
	df_FPKM = subset(df_FPKM,gene %in% df_event$gene)
	df_FPKM = df_FPKM[,!(names(df_FPKM) %in%c("genename","cate","LENGTH"))]
	df_FPKM_melt = melt(df_FPKM)
	colnames(df_FPKM_melt) = c("gene","exp","FPKM")
	print(df_FPKM_melt[1:5,])
	df_event_FPKM = merge(df_event,df_FPKM_melt,by=c("gene","exp"),all.x=T)
	#return(df_event_FPKM)
	dt_event_FPKM = data.table(df_event_FPKM)
	dt_event_summ = dt_event_FPKM[,list(crna_reads=sum(total),FPKM=FPKM[1]),by=c("gene","sample")]
	#return(dt_event_summ)
	return(dt_event_FPKM)
}

cRNA_counts_gene_FPKM_join_other_infos = function(channel,counts_FPKM_table,infotable,event_exp_table,cut){
	sql=sprintf("select a.*,b.*,c.total as TOTAL_READS from %s a join %s b join %s c on a.event=b.event and a.event=c.event where c.total>=%s",counts_FPKM_table,infotable,event_exp_table,cut)
	print(sql)
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
	df = data.frame(transc=c(),in_id=c(),count=c(),len=c(),name=c())
	for(x in 1:length(datanames)){
		dt = datalist[[x]]
		sample = datanames[x]
		result = data.frame(dt[,list(count=sum(family==repeat_type),
                   len=intron_end[1]-intron_start[1]),
             by=c("transc","in_id")])
		print(result[1:10,])
		print(sample)
		result$typename = sample
		df = rbind(df,result[,1:5])} 
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
	blast_infos = cbind(id=1:length(blast_infos[,1]),blast_infos)
	#print(blast_infos[1,])
	#print(blast_infos[1,c('id','event','identity','align_len','mismatches','gap','evalue','bit_score','transc','strand','chr','rc_in1_s','rc_in1_e','rc_in2_s','rc_in2_e')])
	sqlSave(channel,blast_infos[,c('id','event','identity','align_len','mismatches','gap','evalue','bit_score','transc','strand','chr','rc_in1_s','rc_in1_e','rc_in2_s','rc_in2_e')],tablename=out_table,rownames=F,colnames=F)
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

intron_pair_repeat_analysis = function(channel,event_table,repeat_up_table,repeat_down_table,out_table){

	sql1 = sprintf("SELECT a.transc,a.in_id,a.strand AS strand1,b.strand AS strand2,a.family,a.chr,a.start AS s1,a.end AS e1,b.start AS s2,b.end AS e2 
			FROM %s a JOIN %s b ON a.transc=b.transc AND a.in_id=b.in_id AND a.family=b.family AND a.strand!=b.strand WHERE a.end<b.start",repeat_up_table,repeat_up_table)
	sql2 = sprintf("SELECT a.transc,a.in_id,a.strand AS strand1,b.strand AS strand2,a.family,a.chr,a.start AS s1,a.end AS e1,b.start AS s2,b.end AS e2 
			FROM %s a JOIN %s b ON a.transc=b.transc AND a.in_id=b.in_id AND a.family=b.family AND a.strand!=b.strand WHERE a.end<b.start",repeat_down_table,repeat_down_table)
	events = sqlQuery(channel,sprintf("select * from %s",event_table))
	up_self = sqlQuery(channel,sql1)
	down_self = sqlQuery(channel,sql2)

	up_self$l_ex = up_self$in_id+1
	up_self$dist = up_self$s2-up_self$e1
	up_self_info = data.table(up_self[,c("transc","l_ex","dist")])
	up_self_info_summ = up_self_info[,list(updist_min=min(dist), up_k1=sum(dist<1000),up_k2=sum(dist<2000),up_k5=sum(dist<5000)),by=c("transc","l_ex")]
	down_self$r_ex = down_self$in_id
	down_self$dist = down_self$s2-down_self$e1
	down_self_info = data.table(down_self[,c("transc","r_ex","dist")])
	down_self_info_summ = down_self_info[,list(downdist_min=min(dist),down_k1=sum(dist<1000),down_k2=sum(dist<2000),down_k5=sum(dist<5000)),by=c("transc","r_ex")]

	events_1 = merge(events,up_self_info_summ,all.x=T,by=c('transc','l_ex'))
	events_2 = merge(events_1,down_self_info_summ,all.x=T,by=c('transc','r_ex'))

	all_report2 = sqlQuery(channel,sprintf("SELECT a.*,b.start,b.end,b.strand,c.start,c.end,c.strand FROM %s a, %s b,%s c 
	                       WHERE a.`transc`=b.`transc` AND a.`transc`=c.`transc` AND a.`l_ex`=b.`in_id`+1 AND a.`r_ex`=c.`in_id` AND b.`family`=c.`family`",event_table,repeat_up_table,repeat_down_table))

	all_report_sel2 = subset(all_report2,strand.1!=strand.2 & family %in% c("Alu","B2","B4"))
	#all_report_sel_dt = data.table(all_report_sel2[,c(1:10,22,17:19,31:33)])
	#all_report_sel_dt = data.table(all_report_sel2[,c(1:11,23,18:20,32:34)])
	all_report_sel_dt$up_dist=0
	all_report_sel_dt$down_dist=0
	pos = all_report_sel_dt$strand=="+"
	all_report_sel_dt$up_dist[pos]=(all_report_sel_dt$in1_e-all_report_sel_dt$end)[pos]
	all_report_sel_dt$down_dist[pos]=(all_report_sel_dt$start.1-all_report_sel_dt$in2_s)[pos]
	pos = all_report_sel_dt$strand=="-"
	all_report_sel_dt$up_dist[pos]=(all_report_sel_dt$start-all_report_sel_dt$in1_s)[pos]
	all_report_sel_dt$down_dist[pos]=(all_report_sel_dt$in2_e-all_report_sel_dt$end.1)[pos]

	all_report_sel_dt=data.table(all_report_sel_dt)
	all_report_sel_dt_summ = all_report_sel_dt[,list(RC_count=length(l_ex),min_up=min(up_dist),min_down=min(down_dist),min_all=min(up_dist+down_dist)),by=c("event")]
	merged = merge(events_2[,c(4,1:3,5:16)],all_report_sel_dt_summ,all.x=T,by=c("event"))
	sqlSave(channel,merged,out_table,rownames=F,colnames=F)
}


circ_splicing_site_analysis = function(){
	for (x in 2:5){
  	rs[,x] = substr(rs[,x],3,6) 
 	cs[,x] = substr(cs[,x],3,6)}
	rs[(!rs[,2] %in% type5),2] = 'others'
	rs[(!rs[,4] %in% type5),4] = 'others'
	rs[(!rs[,3] %in% type3),3] = 'others'
	rs[(!rs[,5] %in% type3),5] = 'others'
	cs[(!cs[,2] %in% type5),2] = 'others'
	cs[(!cs[,4] %in% type5),4] = 'others'
	cs[(!cs[,3] %in% type3),3] = 'others'
	cs[(!cs[,5] %in% type3),5] = 'others'
	out = data.frame(type3=type3,type5=type5)
	out$rand_up5 = sapply(type5,function(x) sum(rs[,2]==x)/length(rs[,1]))
	out$circ_up5 = sapply(type5,function(x) sum(cs[,2]==x)/length(cs[,1]))
	out$rand_down5 = sapply(type5,function(x) sum(rs[,4]==x)/length(rs[,1]))
	out$circ_down5 = sapply(type5,function(x) sum(cs[,4]==x)/length(cs[,1]))
	out$rand_up3 = sapply(type3,function(x) sum(rs[,3]==x)/length(rs[,1]))
	out$circ_up3 = sapply(type3,function(x) sum(cs[,3]==x)/length(cs[,1]))
	out$rand_down3 = sapply(type3,function(x) sum(rs[,5]==x)/length(rs[,1]))
	out$circ_down3 = sapply(type3,function(x) sum(cs[,5]==x)/length(cs[,1]))
}

cRNA_flanking_intron_repeats = function(channel,database,inevent,inrepeat,in1_repeats,in2_repeats){
	#sqlQuery(channel,paste("CREATE TABLE ",database,"temp_r1 AS SELECT transc,l_ex-1 AS in_id,chr,strand,in1_s AS `start`,in1_e AS `end` FROM ",database,inevent,' GROUP by transc,in_id',sep=""))
	sqlQuery(channel,paste("CREATE TABLE ",database,"temp_r2 AS SELECT transc,r_ex AS in_id,chr,strand,in2_s AS `start`,in2_e AS `end` FROM ",database,inevent,' GROUP by transc,in_id',sep=""))	
	#sqlQuery(channel,paste('ALTER TABLE ',database,'temp_r1 ADD KEY (`chr`,`start`,`end`)',sep=""))
	sqlQuery(channel,paste('ALTER TABLE ',database,'temp_r2 ADD KEY (`chr`,`start`,`end`)',sep=""))
	sql1 = sprintf("create table %s%s as SELECT a.transc,a.`in_id`,a.`strand` AS intron_strand,a.`start` AS intron_start,a.`end` AS intron_end,b.* FROM %stemp_r1 a JOIN %s b ON a.`chr`=b.`chr` AND b.`start`>a.`start` AND b.`end`<a.`end`",
		database,in1_repeats,database,inrepeat)
	sql2 = sprintf("create table %s%s as SELECT a.transc,a.`in_id`,a.`strand` AS intron_strand,a.`start` AS intron_start,a.`end` AS intron_end,b.* FROM %stemp_r2 a JOIN %s b ON a.`chr`=b.`chr` AND b.`start`>a.`start` AND b.`end`<a.`end`",
		database,in2_repeats,database,inrepeat)
	print(sql1)
	print(sql2)
	#sqlQuery(channel,sql1)
	sqlQuery(channel,sql2)
	#sqlQuery(channel,sprintf("alter table %s%s add key (`transc`,`in_id`,`family`)",database,in1_repeats))
	sqlQuery(channel,sprintf("alter table %s%s add key (`transc`,`in_id`,`family`)",database,in2_repeats))
}


