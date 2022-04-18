
setwd("~/Desktop/Covid_Atlas/Pseudotime/Integrated_lung_epithelial/")
###lung epithelial is the merged dataset by donor
dge=list()
for (i in 1:length(unique(lung_epithelial$X10X_run))){
dge[[i]] <- subset(lung_epithelial,cells=colnames(lung_epithelial)[lung_epithelial$X10X_run==unique(lung_epithelial$X10X_run)[i]])
# dge[[i]] <- subset(lung_epithelial,cells=colnames(lung_epithelial)[lung_epithelial$donor==unique(lung_epithelial$donor)[i]])
}

for (i in 1:length(dge)) {
dge[[i]] <- NormalizeData(dge[[i]], verbose = FALSE)
dge[[i]] <- FindVariableFeatures(dge[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}
dge.anchors <- FindIntegrationAnchors(object.list = dge, dims = 1:30)
dge <- IntegrateData(anchorset = dge.anchors, dims = 1:30,k.weight = 30)

###normalize RNA assay for gene expression comparison

DefaultAssay(dge) <- 'RNA'
genes.use = Seurat:::CheckFeatures(object.name = dge,data  = dge[['RNA']]@data, features = rownames(dge[['RNA']]@data),verbose = F)
dge =  ScaleData(dge, features = genes.use,vars.to.regress = c('nCount_RNA'))


DefaultAssay(dge) <- "integrated"

# Run the standard workflow for visualization and clustering
dge <- ScaleData(dge, verbose = FALSE)
dge <- RunPCA(dge, npcs = 100, verbose = FALSE)
# t-SNE and Clustering
dge <- RunUMAP(dge, reduction = "pca", dims = 1:100)
dge <- FindNeighbors(dge, reduction = "pca", dims = 1:100)
dge <- FindClusters(dge, resolution = 0.5)

# Visualization

DimPlot(dge, reduction = "umap", group.by = "X10X_run",label = F)
ggsave(file="UMAP_integrated_bydonor.pdf",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "umap", group.by = "cell_type_annotation",label = T)
ggsave(file="UMAP_integrated_byID.pdf",width = 20,height = 20,units = "cm")
DimPlot(dge, reduction = "umap", label = T)
ggsave(file="UMAP_integrated_bycluster.pdf",width = 20,height = 20,units = "cm")

DimPlot(dge, reduction = "umap", label = T,group.by = "ID")
ggsave(file="UMAP_integrated_byID.pdf",width = 20,height = 20,units = "cm")

DimPlot(dge, reduction = "umap", label = F,group.by = "X10X_run")
ggsave(file="UMAP_integrated_bydonor.pdf",width = 20,height = 20,units = "cm")

#####generate signature
dge<-SetIdent(dge,value = as.vector(dge@active.ident))
for (i in 1:length(unique(dge@active.ident))) {
  markers_temp <- FindMarkers(dge, ident.1 = unique(dge@active.ident)[i], ident.2 = NULL, only.pos = TRUE)
  markers_temp$gene<-as.vector(rownames(markers_temp))
  markers_temp<-markers_temp[order(markers_temp$avg_log2FC,decreasing = T),]
  sig_cluster<-markers_temp$gene[markers_temp$p_val_adj<0.05]
  sig_cluster<-sig_cluster[1:min(length(sig_cluster),50)]
  write.table(sig_cluster,paste0("",unique(dge@active.ident)[i],"_Signature.txt"),sep = "\t",quote = F,col.names = F,row.names = F)
}


###Make heatmap with dendrogram
library(tidyr)
library(textshape)
pbmc.markers <- FindAllMarkers(dge,assay = "RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
genelist<-as.vector(top20$gene)
dge_temp<-subset(dge,downsample=500)
z <- DoHeatmap(dge_temp, features = genelist) + NoLegend()
y <- z$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()

my_sample_col <- data.frame((dge$cell_type_annotation))
colnames(my_sample_col)="Annotation"
my_gene_col <-as.data.frame(c(rep("AT1_Marker",20),rep("AT2_Marker",20),rep("Basal_Marker",20),rep("Transitional_Marker",20)))
rownames(my_gene_col)=t(make.unique(genelist))
colnames(my_gene_col)="Markers"
pheatmap::pheatmap(w,annotation_row = my_gene_col,annotation_col = my_sample_col,show_colnames = F,clustering_method="ward.D2")




