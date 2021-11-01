library(DropletUtils)
library(Matrix)
library(Seurat)
library(dplyr)

# read from .mtx (better)
counts = readMM('sparse_matrix.mtx')
# obs data.frame
obs = read.csv('obs.csv', header = T)
# var data.frame
vars = read.csv('var.csv', row.names = 1)
# obsm UMAP coords
obsm = read.csv('obsm.csv', header = T)

#check that vars is not empty
if(dim(vars)[2]==0){
    vars$gene_symbol = row.names(vars)
    #vars$gene_ids = row.names(vars)
    #vars$feature_types = 'Gene Expression'
}

# for CTA datasets, we replace the name of the virus transcripts
# fix the virus name
#vars$gene_symbol[grep('SARS',vars$gene_ids)] <- vars$gene_ids[grep('SARS',vars$gene_ids)]
#fix the exmpy first element in vars (some bug from scanpy after concatenation)
if(vars$gene_symbol[1]=="") 
    vars$gene_symbol[1] = "101"

colnames(counts) <- vars$gene_symbol
row.names(counts) <- obs[,1]

row.names(obs) <- obs[,1]
row.names(obsm) <- obs[,1]

# rename as cell id
#obs %>% rename(cell_id = X) -> obs

scvi_umap = obsm[, c('X_umap1', 'X_umap2')]
row.names(scvi_umap) <- row.names(obs)

# Make Seurat object
seurat_obj = CreateSeuratObject(counts = t(counts), meta.data = obs)


seurat_obj[['scvi_umap']] =CreateDimReducObject(embeddings = scvi_umap %>% as.matrix() , key="SCVI_", assay =DefaultAssay(seurat_obj))

# save to .rds
save(seurat_obj,file = 'seurat_obj.rds')
