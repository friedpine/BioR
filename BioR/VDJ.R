get_ranges_mean_depth =  function(channel,table_ranges,table_depth){
	sql = paste('SELECT a.*,b.sample,SUM(b.depth)/(`end`-`start`) AS mean_depth FROM ',table_ranges,
		' a JOIN ',table_depth,
		' b ON b.pos>=a.`start` AND b.pos<=a.`end` GROUP BY b.sample,a.gene order by a.`start`')
	print(sql)
	range_info = sqlQuery(channel,sql)
	genes = as.character(range_info$gene[!duplicated(range_info$gene)])
	samples = as.character(range_info$sample[!duplicated(range_info$sample)])
	out = matrix(0,ncol=length(samples),nrow=length(genes))
	out_len = length(range_info$gene)
	colnames(out) = samples
	rownames(out) = genes
	for (x in 1:out_len){
		out[range_info[x,'gene'],range_info[x,'sample']]=range_info[x,'mean_depth']}
	
	return(t(out))
}

get_range_IGH_overlaption = function(channel,table_range,table_info){
	sql = paste('SELECT * FROM ',table_range,' a JOIN ',table_info,' b ON ',   
						'(a.range_up>b.start-50 AND a.range_up<b.end+50) OR
						(a.range_down>b.start-50 AND a.range_down<b.end+50) OR
						(b.start>=a.range_up AND b.end<=a.range_down)')
	range_info = sqlQuery(channel,sql)
	return(range_info)
}

get_junction_IGH_overlaption = function(channel,table_range,table_info){
	sql = paste('SELECT * FROM ',table_range,' a JOIN ',table_info,' b ON ',   
						'(a.range_up>b.start-50 AND a.range_up<b.end+50) OR
						(a.range_down>b.start-50 AND a.range_down<b.end+50) OR
						(b.start>=a.range_up AND b.end<=a.range_down)')
	range_info = sqlQuery(channel,sql)
	return(range_info)
}



