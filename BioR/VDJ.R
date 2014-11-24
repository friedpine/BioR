get_ranges_mean_depth =  function(channel,table_ranges,table_depth){
	sql = paste('SELECT a.*,AVG(b.depth) AS mean_depth FROM ',table_ranges,
		' a JOIN ',table_depth,
		' b ON a.sample=b.sample AND b.pos>=a.range_up AND b.pos<=a.range_down GROUP BY a.id')
	range_info = sqlQuery(channel,sql)
	return(range_info)
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

