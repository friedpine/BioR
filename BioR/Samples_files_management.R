library(RODBC)

#INPUT1: files lists and informations file (path & size)
#INPUT2: files lists RAW_sample names and official names maybe other infos(rawname & sample)

fastq_file_lists = function(samples,files,types){
  files = files[,c(5,9)]
  colnames(files) = c("size","path")
  files$rawname = ldply(strsplit(as.character(files$path),"/"),function(x) x[length(x)-1])[,1]
  files$type = types
  merged = merge(files,samples,by="rawname")
  df_up = data.frame(sample=merged$sample,type=merged$type,
                    path=merged$path,method="A",total_reads=10000000)
  return(df_up)
}

#channel <- odbcConnect("IBM", uid="root", pwd="123456")

#FUCTION: Save the lists of samples to the database!
#Requirement: 
fastq_to_database = function(channel,database,table,data){
table = paste(database,".",table,sep="")
sql = paste("CREATE TABLE ",table," (
  `sample` varchar(255) DEFAULT NULL,
  `type` varchar(255) DEFAULT NULL,
  `path` varchar(1000) DEFAULT NULL,
  `size` double DEFAULT NULL,
  `method` text,
  UNIQUE KEY `u` (`sample`,`type`))",sep="")
show_tables = sqlTables(channel,catalog=database)
if (!(table %in% show_tables$TABLE_NAME)){sqlQuery(channel,sql)}
df_upload = data[,c("sample","type","path","size")]
df_upload$method = "NA"
print(df_upload[1,])
sqlSave(channel,dat=df_upload,tablename=table,append=T,rownames=F,colnames=F)
#return(df_upload)
} 




