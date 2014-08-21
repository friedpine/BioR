#INPUT1: files lists and informations file (path & size)
#INPUT2: files lists RAW_sample names and official names maybe other infos(rawname & sample)

fastq_file_lists = function(samples,files,types){
  files$rawname = ldply(strsplit(as.character(files$path),"/"),function(x) x[length(x)-1])
  files$type = types
  merged = merge(files,samples,by="rawname")
  df_up = data.frame(sample=merged$sample,type=merged$type,
                    path=merged$path,method="A",total_reads=1000000)
  return(df_up)
}
