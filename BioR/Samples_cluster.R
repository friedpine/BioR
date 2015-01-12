library(FactoMineR)

#USING_PCA_METHOD
sample_PCA = function(data,low_cut,low_log2,expressed_cells,dimensions,genes_numbers,cell_types){
  data.exp = data[rowSums(data>=low_cut)>=expressed_cells,]
  data.exp.log2 = log2(data.exp+low_log2)
  PCA1 = PCA(t(data.exp.log2))
  dimension.PCA.allgenes<-dimdesc(PCA1, axes=c(1,2,3,4))
  df_PCA_score = data.frame(PC1=PCA1$ind$coord[,1],PC2=PCA1$ind$coord[,2],
                  PC3=PCA1$ind$coord[,3],PC4=PCA1$ind$coord[,4],
                type=cell_types)
  #get_HIGH_LOADING_VALUE_genes
  dim4<-as.data.frame(dimension.PCA.allgenes[[4]])
  dim3<-as.data.frame(dimension.PCA.allgenes[[3]])
  dim2<-as.data.frame(dimension.PCA.allgenes[[2]])
  dim1<-as.data.frame(dimension.PCA.allgenes[[1]])
  high_loading_genes = unique(c(rownames(dim1[order(abs(dim1$quanti.correlation),decreasing=T)[1:genes_numbers],]),
                                      rownames(dim2[order(abs(dim2$quanti.correlation),decreasing=T)[1:genes_numbers],]),
                                      rownames(dim3[order(abs(dim3$quanti.correlation),decreasing=T)[1:genes_numbers],]),
                                      rownames(dim4[order(abs(dim4$quanti.correlation),decreasing=T)[1:genes_numbers],])))  
  return(list(genes=data.exp.log2[high_loading_genes,],score=df_PCA_score))
  }

#CLUSTERING_PLOTTING
sample_cluster_plot_infos = function(PCA_scores,clustering_data,geneID2name,filename){
  library(ggplot2)
  library(gplots)
  pdf(filename)

  print(ggplot(PCA_scores,aes(PC1,PC2,color=type))+
    geom_point(size=5))
    print(ggplot(PCA_scores,aes(PC1,PC3,color=type))+
    geom_point(size=5))
    print(ggplot(PCA_scores,aes(PC1,PC4,color=type))+
    geom_point(size=5))

  heat_data = as.matrix(clustering_data)
  dd_gene <- as.dendrogram(hclust(as.dist((1 - cor(t(heat_data)))/2),))
  dd_gene = reorder(dd_gene,wts=rowMeans(heat_data))
  order_gene = sapply(cut(dd_gene,0)[[2]],function(x) attr(x,'label'))
  dd_sample <- as.dendrogram(hclust(dist(t(heat_data),method = "euclidean")))
  dd_sample = reorder(dd_sample,wts=colMeans(heat_data))
  order_sample = sapply(cut(dd_sample,0)[[2]],function(x) attr(x,'label'))
  heatmap.2(heat_data,dendrogram="both",
          scale="row",density.info='none',Rowv=dd_gene,Colv=dd_sample,col=bluered(100),trace='none')
  dev.off()

  rownames(geneID2name) = geneID2name$gene
  outtable = as.data.frame(heat_data[rev(order_gene),order_sample])
  outtable$gene = rownames(outtable)
  outtable$genename = geneID2name[outtable$gene,'genename']
  write.csv(outtable,file=gsub('pdf','csv',filename))
}

#CLUSTERING_PLOTTING
sample_cluster_plot = function(PCA_scores,clustering_data,filename){
  library(ggplot2)
  library(gplots)
  pdf(filename)

  print(ggplot(PCA_scores,aes(PC1,PC2,color=type))+
    geom_point(size=5))
    print(ggplot(PCA_scores,aes(PC1,PC3,color=type))+
    geom_point(size=5))
    print(ggplot(PCA_scores,aes(PC1,PC4,color=type))+
    geom_point(size=5))

  heat_data = as.matrix(clustering_data)
  dd_gene <- as.dendrogram(hclust(as.dist((1 - cor(t(heat_data)))/2),))
  dd_gene = reorder(dd_gene,wts=rowMeans(heat_data))
  order_gene = sapply(cut(dd_gene,0)[[2]],function(x) attr(x,'label'))
  dd_sample <- as.dendrogram(hclust(dist(t(heat_data),method = "euclidean")))
  dd_sample = reorder(dd_sample,wts=colMeans(heat_data))
  order_sample = sapply(cut(dd_sample,0)[[2]],function(x) attr(x,'label'))
  heatmap.2(heat_data,dendrogram="both",
          scale="row",density.info='none',Rowv=dd_gene,Colv=dd_sample,col=bluered(100),trace='none')
  dev.off()
}
