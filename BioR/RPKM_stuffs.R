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

XPKM_corr_scatter_matrix = function(data,samples1,samples2,low_cut,xlimit,plotfile){
  samples = colnames(data)
  counts = length(samples)
  cormat = matrix(0,nrow=counts,ncol=counts)
  colnames(cormat) = samples
  rownames(cormat) = samples
  pdf(plotfile)
  for (x in samples1){
    for (y in samples2){
      data_exp = data[,c(x,y)]
      data_nonzero = data_exp[rowSums(data_exp>low_cut)>0,]
      cormat[x,y] = round(cor(data_nonzero[,1],data_nonzero[,2]),2)
      cormat[y,x] = cormat[x,y]
      title_name = paste(samples[x],samples[y],sep="_")
      print(ggplot(data.frame(x=data_nonzero[,1],y=data_nonzero[,2]),aes(x,y))+
        geom_point(size=1)+
        xlim(xlimit[1],xlimit[2])+ylim(xlimit[1],xlimit[2])+
        labs(title = title_name)
        )
    }
    cormat[x,x] = 1
  }
  cormat[counts,counts] = 1
  dev.off()
  return(cormat)
}

corr_coef_matrix = function(data,low_cut,mtd){
	samples = colnames(data)
	counts = length(samples)
	cormat = matrix(0,nrow=counts,ncol=counts)
	colnames(cormat) = samples
	rownames(cormat) = samples
	for (x in seq(1,counts-1)){
		print(x)
		for (y in seq(x+1,counts)){
			data_exp = data[,c(x,y)]
			data_nonzero = data_exp[rowSums(data_exp>low_cut)>1,]
			cormat[x,y] = round(cor(data_nonzero[,1],data_nonzero[,2],method = mtd),2)
			cormat[y,x] = cormat[x,y]
		}
		cormat[x,x] = 1
	}
	cormat[counts,counts] = 1
	return(cormat)
}

#USING_LINEAR_REGRESSION_MODEL!
esitimate_mole_count_with_ERCC = function(channel,exp_table,ERCC_detail,samples_names,sample_ERCC_counts,filename){
	sqlRPKM = paste("select genename,",toString(samples_names,sep=",")," from ",exp_table,sep="")
	sqlERCC = paste("select * from ",ERCC_detail,sep="")
	print(sqlRPKM)
	data_exp = sqlQuery(channel,sqlRPKM)
	data_ercc = sqlQuery(channel,sqlERCC)
	print(data_exp[1,])
	data_exp_ercc = subset(data_exp,genename %in% data_ercc[,1])
	data_exp_mRNA = subset(data_exp,!(genename %in% data_ercc[,1]))
	out = array()
	ERCC_reads = array()
	print(data_exp_ercc[1,])
	pdf(filename)
	for (x in 1:length(samples_names)){
		sample = samples_names[x]
		print(sample)
		ERCC_total = sample_ERCC_counts[x]
		scale_coeff = ERCC_total/(70000000)
		exp_ercc = data_exp_ercc[,c(1,x+1)]
		exp_mRNA = data_exp_mRNA[,x+1]
		count_ercc = data_ercc
		count_ercc$count = ERCC_total*count_ercc$count
		df_ercc_count_exp = merge(exp_ercc,count_ercc,by="genename")
		colnames(df_ercc_count_exp) = c("gene","exp","count")
		#out[[x]]=df_ercc_count_exp
		df_ercc_count_exp$log10count = log10(df_ercc_count_exp$count+0.01)
		df_ercc_count_exp$log10exp = log10(df_ercc_count_exp$exp+0.01)
		data=subset(df_ercc_count_exp,log10exp>2+log10(scale_coeff)&log10count>4+log10(scale_coeff))
		print(df_ercc_count_exp[1:10,])
		print(data[1,])
		model = lm(log10count~log10exp,data=data)

		#Predicted the whole ERCC genes including the filtered!
		new_ercc = data.frame(log10exp=df_ercc_count_exp$log10exp)
		df_ercc_count_exp$predict = 10**(predict(model, newdata = new_ercc))
		

		new_mRNA = data.frame(exp=exp_mRNA,log10exp=log10(exp_mRNA+0.01))
		new_mRNA$predict = 10**(predict(model, newdata = new_mRNA))
		print(colSums(new_mRNA))
		out[x] = colSums(new_mRNA)[3]
		print(ggplot(df_ercc_count_exp,aes(log10exp,log10count))+
		  geom_point()+
		  geom_abline(intercept=model$coefficients[1],slope=model$coefficients[2])+
		  labs(title = sample,x="log10(RPKM)",y="log10(Molecular_count)"))		
	}
	dev.off()
	return(data.frame(sample=samples_names,pred_count=out))
}

#USING_LINEAR_REGRESSION_MODEL_TO_Predict_circularRNA_molecular_Counts!!!
esitimate_circ_count_with_ERCC = function(channel,data_exp_mRNA,exp_table,ERCC_detail,samples_names,sample_ERCC_counts,filename){
	sqlERCC = paste("select * from ",ERCC_detail,sep="")
	sqlRPKM = paste("select genename,",toString(samples_names,sep=",")," from ",exp_table,sep="")
	data_ercc = sqlQuery(channel,sqlERCC)
	data_exp = sqlQuery(channel,sqlRPKM)	
	data_exp_ercc = subset(data_exp,genename %in% data_ercc[,1])

	out = array()
	ERCC_reads = array()
	for (x in 1:length(samples_names)){
		sample = samples_names[x]
		print(sample)

		#Prepare circular prediction data!
		exp_mRNA = data_exp_mRNA[,sample]
		new_mRNA = data.frame(exp=exp_mRNA,log10exp=log10(exp_mRNA+0.01))

		#ERCC name/Expression level/molecular Counts!!
		ERCC_total = sample_ERCC_counts[x]
		scale_coeff = ERCC_total/(70000000)
		exp_ercc = data_exp_ercc[,c(1,x+1)]
		
		count_ercc = data_ercc
		count_ercc$count = ERCC_total*count_ercc$count
		df_ercc_count_exp = merge(exp_ercc,count_ercc,by="genename")
		colnames(df_ercc_count_exp) = c("gene","exp","count")
		
		#ERCC expression/counts are log10!!!
		df_ercc_count_exp$log10count = log10(df_ercc_count_exp$count+0.01)
		df_ercc_count_exp$log10exp = log10(df_ercc_count_exp$exp+0.01)
		
		#Then are subsetted! And linera regression!
		data=subset(df_ercc_count_exp,log10exp>2+log10(scale_coeff)&log10count>4+log10(scale_coeff))
		print(df_ercc_count_exp[1:10,])
		print(data[1,])
		model = lm(log10count~log10exp,data=data)
		
		#Predict the mRNA genes!! 
		new_mRNA$predict = 10**(predict(model, newdata = new_mRNA))
		print(colSums(new_mRNA))
		out[x] = colSums(new_mRNA)[3]
	}
	dev.off()
	return(data.frame(sample=samples_names,pred_circ_count=out))
}