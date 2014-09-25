library(ggplot2)
scatter_corr = function(channel,sample_table,table,s1,s2,col_name,low_cut,filename){
	s1 = sqlQuery(channel,paste("select `exp` from ",sample_table," where sample= '",s1,"'",sep=""))
	s2 = sqlQuery(channel,paste("select `exp` from ",sample_table," where sample= '",s2,"'",sep=""))
	data_exp = sqlQuery(channel,paste("select ",s1[1,1],",",s2[1,1]," from ",table,sep=""))
	data_nonzero = data_exp[rowSums(data_exp>low_cut)>0,]
	colnames(data_nonzero) = c("v1","v2")
	corr_coef = round(cor(data_nonzero[,1],data_nonzero[,2]),2)
	print(corr_coef)
	pdf(filename)
	p1 = ggplot(data_nonzero,aes(v1,v2))+
			geom_point(size=1)+
			annotate("text",label=paste("r","=",corr_coef),x=1000,y=8500,size=9)+
			scale_x_log10(limits=c(1, 10000),breaks=c(1,10,100,1000,5000))+
			scale_y_log10(limits=c(1, 10000),breaks=c(1,10,100,1000,5000))+
			xlab(col_name[1])+
			ylab(col_name[2])
	print(p1)
	dev.off()
}
