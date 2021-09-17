# Differential expression using MAST 
# Input: Annotated Seurat object 
#  



library(dplyr)
library(zinbwave)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(MAST)
library(doParallel)
library(BiocParallel)


# Sample ID 
tissue = 'liver'
# Cell type for DGE analysis 
cell_type = 'hepatocyte' 
# filter lowly expressed genes
min_counts = 2
min_cells = 5
# How many genes to consider for DGE (top n based on high variance)
n_top_genes = 0
# output 
output_dir='/mnt/ibm_lg/covid_tissue_atlas/results/DGE/'

#### 1. Prepare data 
# Assumes general data structure exported from scanpy 
# read from .mtx (better)
counts = readMM('sparse_matrix.mtx')
# obs data.frame
obs = read.csv('obs.csv', header = T)
# var data.frame
vars = read.csv('var.csv')
# obsm UMAP coords
obsm = read.csv('obsm.csv', header = T)

# Covariates for DGE 
covariate = ' + lab'
covariate = ''


# for CTA datasets, we replace the name of the virus transcripts
# fix the virus name
if(grep(names(vars), pattern = 'gene_id')){
  vars$gene_symbol[grep('SARS',vars$gene_ids)] <- vars$gene_ids[grep('SARS',vars$gene_ids)]
}

colnames(counts) <- vars$gene_symbol
row.names(counts) <- obs[,1]

row.names(obs) <- obs[,1]
row.names(obsm) <- obs[,1]

# Create SCE object 
# for some reason is not default in readMM
c <- as(counts, 'dgCMatrix') %>% as.matrix
sce <- SingleCellExperiment(list(counts=t(c)),
                            colData=obs,
                            rowData=vars,
                            metadata=list(tissue= tissue)
)


# 2. Set up metadata for differntial gene expression
core = sce[,colData(sce)$annotations_v1 %in% c(cell_type )]
# relaxed gene filter: 
core = core[rowSums(assay(core)) > 0, ]
# stronger filter: only genes expressed 
# ideally we want to use highly variable genes. 
# but we don't want to compute DGE from genes expressed in a few cells. 
filter <- rowSums(assay(core)>min_counts)>min_cells
core <- core[filter,]

# further filter by high variance 
assay(core) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(core)
# top variable genes -- quick and dirty tbh
vars <- sort(vars, decreasing = TRUE)
if(n_top_genes){
    core <- core[names(vars)[1:n_top_genes], ]
}

# setup for model 
#colData(core)$disease_status = factor(colData(core)$disease_status)
 

NCORES = 25 
registerDoParallel(NCORES)
register(DoparParam())



# Export count matrix
counts = assay(core)

# MAST requires log data
tpm <- counts*1e6/colSums(counts)
tpm <- log2(tpm+1)

core_log = core 
assay(core_log) <- tpm

# MAST sca object structure
sca <- SceToSingleCellAssay(core_log, class = "SingleCellAssay")


zlmCond <- zlm(formula = as.formula(paste0("~disease_status", covariate)), sca=sca)

summaryCond <- summary(zlmCond, doLRT="disease_statusCov19")

summaryDt <- summaryCond$datatable
dt1 = summaryDt[contrast=="disease_statusCov19" & component=="H", .(primerid, `Pr(>Chisq)`)]
dt2 = summaryDt[contrast=="disease_statusCov19" & component=="logFC", .(primerid, coef, z)]
de_res = merge(dt1, dt2, by="primerid")
colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "BH")

# convert to iDEA 
mast_idea <- function(res = data.frame() , use_mast_zscore = T){
    
    if(use_mast_zscore){
        zscore <- res$age.logFC_z
    }else{   
        pvalue <- res$age.H_fdr
        zscore <- qnorm(pvalue/2.0, lower.tail = FALSE)
    }
    beta <- res$age.logFC
    se_beta <- abs(beta/zscore)
    beta_var <- se_beta^2
    summary_deseq  = data.frame(beta = beta, beta_var = beta_var)
    rownames(summary_deseq) <- res$gene
    return(summary_deseq)
}

# get summary data frame 
summary_mast <- mast_idea(res= de_res)

# write output 
write.csv(summary_mast, file =paste0(output_dir, "DGE_", tissue,"_",cell_type,".csv"))