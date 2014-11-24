library(data.table)
library(ggplot2)
library(reshape2)
library(gplots)
library(plyr)
library(RODBC)

circRNA_events_exp_infos = function(channel,intron_table,exonmaps_tables,out_events,out_eventsinfo,out_exp){
	all_events = data.frame(Sample=character(),
                 event=character(), 
                 transc=character(), 
                 l_ex=integer(),
                 r_ex=integer(),
                 stringsAsFactors=FALSE) 
	for (x in exonmaps_tables){
		temp = sqlQuery(channel,paste("select sample,event,transc,l_ex,r_ex from ",x," where t_dis=0 GROUP BY sample,`event`,l_dis"))
		if(!is.data.frame(temp)){print("FALSE QUERY!");return(0);}
		all_events = rbind(all_events,temp)
	}
	dt = data.table(all_events)
	circevents = dt[,list(transc=transc[1],l_ex=l_ex[1],r_ex=r_ex[1],total=length(sample)),by=c("event")]
	circevents_exp = dt[,list(total=length(l_ex),transc=transc[1]),by=c("event","sample")]
	sqlSave(channel,circevents,tablename=out_events,rownames=F,colnames=F)
	sqlSave(channel,circevents_exp,tablename=out_exp,rownames=F,colnames=F)
	sqlQuery(channel,paste("alter table ",out_events," add index tsc1 (transc,l_ex)"))
	sqlQuery(channel,paste("alter table ",out_events," add index tsc2 (transc,r_ex)"))
	info_l = sqlQuery(channel,paste("select * from ",out_events," a left join ",intron_table," b on a.transc=b.transc and a.l_ex=b.ex_order+1"))
	info_r = sqlQuery(channel,paste("select * from ",out_events," a left join ",intron_table," b on a.transc=b.transc and a.r_ex=b.ex_order"))
	info_out = cbind(info_l[,c(1,2,3,4,5,8,9,10,11)],info_r[,c(10,11)])
	colnames(info_out) = c("event","transc","l_ex","r_ex","exp","chr","strand",'in1_s','in1_e','in2_s','in2_e')
	sqlSave(channel,info_out,tablename=out_eventsinfo,rownames=F,colnames=F)
}

circ_info_from_transc_pos = function(channel,eventtable,intron_table,outtable){
	sql = paste("create table",outtable," as SELECT a.event,a.transc,a.l_ex,a.`r_ex`,b.`chr`,b.`strand`,b.`left_pos` AS in1_s,b.`right_pos` AS in1_e,c.`left_pos` AS in2_s,c.`right_pos` AS in2_e FROM",
	 eventtable,"a JOIN",intron_table,"b JOIN",intron_table,"c ON a.transc=b.transc AND a.transc=c.transc AND a.`l_ex`=b.`ex_order`+1 AND  a.`r_ex`=c.`ex_order`",sep=" ")
	sqlQuery(channel,sql)
}


