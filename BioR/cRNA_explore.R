library(data.table)
library(ggplot2)

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
