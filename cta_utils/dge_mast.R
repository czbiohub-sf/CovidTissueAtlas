# Differential expression using MAST 
# Input: Annotated Seurat object 
#  


suppressPackageStartupMessages({
    library(dplyr)
    library(zinbwave)
    library(SingleCellExperiment)
    library(Matrix)
    library(Seurat)
    library(MAST)
    library(doParallel)
    library(BiocParallel)

})

DGE_DIR = '/mnt/ibm_lg/covid_tissue_atlas/results/DGE/MAST/'

FCTHRESHOLD = log2(1.5)
# output 
output_dir=

# raw input files 
NCORES = 50
#registerDoParallel(NCORES)
#register(DoparParam())

options(mc.cores = NCORES)


# This function runs DGE on all organs 
# writes a .csv with all DE genes for each organ 
# output dir default: '/mnt/ibm_lg/covid_tissue_atlas/results/DGE/covid_control_celltypes/'
# use data_id to provide specific suffix for each run and prevent overwritting 
dgeCTA <- function(regr_ngenes = T, out_id = ''){
    tissues = c('liver', 'kidney', 'lung',  'heart', 'prostate')

    for(t in 1:length(tissues)){
        dgeOrgan(tissue =tissues[t], ann_col = 'cell_type_annotation', 
                data_id = out_id, 
                regress_ngenes = regr_ngenes)
        print(paste0('DONE tissue ', tissues[t]))
        gc() 
    }
}

# Run all cell types for a given organ 
dgeOrgan <- function(tissue = "liver",  # tissue to process
                    ann_col = "cell_type_annotation", # which col in meta.data has the cell type label
                    output_dir = '/mnt/ibm_lg/covid_tissue_atlas/results/DGE/MAST/', # write DGE .csvs here 
                    input_dir_prefix= '/mnt/ibm_lg/covid_tissue_atlas/data/tissue_objects/all_tissues/', 
                    data_id = "", # ID append for output files  , 
                    regress_ngenes = T
                    ){
    # 1. read from organ directory 
    input_dir= paste0(input_dir_prefix, tissue, "/r_files/")
    sce <- load_sce_object(input_dir, organ = tissue )
    # 2. list all cell types for this organ 
    colData(sce)[ann_col][,1] %>% unique() -> cell_type_list

    all_dge = data.frame() 

    # Regress number of detected genes 
    if(regress_ngenes){
        cov_var = " + cgeneson"
    }else{ 
        cov_var = ""
    }
            
    # 3. DGE for each cell type 
    for(c in 1:length(cell_type_list)){

        sce_sub <- splitAndProcess(sce, ann_labels = ann_col, cell_type = cell_type_list[c], min_counts = 1, min_cells = 5)
        
        n_factors <- colData(sce_sub)['disease_status'][,1] %>% table %>% length
        # run differential expression ONLY if there are 2 factors for this cell type
        if(n_factors>1){
            # NOTE: by default we correct for number of detected genes 

            fc_results <- dge_MAST(sce_sub, covariate = cov_var) 
            # 4. Write files for each cell type
            # collapse white spaces 
            #celltype_file = paste(str_split(cell_type_list[c], " ")[[1]], collapse ="_")
            # remove special characters 
            #celltype_file =  str_replace(celltype_file,'/', '_')
            fc_results$cell_type  = cell_type_list[c]

            all_dge = rbind(all_dge, fc_results)
        }
    }
    write.csv(all_dge, file = paste0(output_dir, tissue,"_DGE_allcelltypes_" , data_id, ".csv") , quote = F, row.names = T)
}

# Reads exported files from scanpy and creates a SCE object 
# counts are assumed to be raw
# Since QC already took place in scanpy, here we only focus on creating the object 
load_sce_object <- function( file_path = "",
                             organ = ""){
    #### 1. Prepare data 
    # Assumes general data structure exported from scanpy 
    # read from .mtx (better)
    counts = readMM(paste0(file_path,'sparse_matrix.mtx'))
    # obs data.frame
    obs = read.csv( paste0(file_path, 'obs.csv'), header = T)
    # var data.frame
    vars = read.csv( paste0(file_path, 'var.csv'), row.names = 1)
    # obsm UMAP coords
    obsm = read.csv( paste0(file_path,'obsm.csv'), header = T)

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

    # Create SCE object 
    # for some reason is not default in readMM
    c <- as(counts, 'dgCMatrix') %>% as.matrix
    sce <- SingleCellExperiment(list(counts=t(c)),
                                colData=obs,
                                rowData=vars,
                                metadata=list(tissue= organ)
    )

    return(sce)

}

# From the organ object, subset for a specific cell type 
# Performs normalization and filtering of lowly expressed genes 
splitAndProcess <- function(sce =c() , ann_labels = 'cell_type_annotation', 
                            cell_type = 'hepatocyte',
                            min_counts = 2, # filter lowly expressed genes
                            min_cells = 5, 
                            n_top_genes = 0 # default, consider all genes in the object for DGE 
                            ){

    # 2. Set up metadata for differntial gene expression
    
    cell_idx <- colData(sce)[ann_labels] %in% c(cell_type) %>% unlist()
    core <- sce[,cell_idx]

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

    # Number of genes expressed (to use as covariate)
    cdr2 <-colSums(assay(core)>0)
    colData(core)$cgeneson <- scale(cdr2)

    # top variable genes -- quick and dirty tbh
    vars <- sort(vars, decreasing = TRUE)
    if(n_top_genes){
        core <- core[names(vars)[1:n_top_genes], ]
    }

    return(core)
}
# setup for model 
#colData(core)$disease_status = factor(colData(core)$disease_status)
 

# # # # # #
 # # # # # #
# MAST 
# for 8k genes this function takes around 10 min on leela 
# returns data.frame of sig DE gens at 0.05 p-value FDR adjusted, and min log2FC 1.5
dge_MAST <- function(core = c(), covariate = "" ){
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
    # dt1 = summaryDt[contrast=="disease_statusCov19" & component=="H", .(primerid, `Pr(>Chisq)`)]
    # dt2 = summaryDt[contrast=="disease_statusCov19" & component=="logFC", .(primerid, coef, z)]
    # de_res = merge(dt1, dt2, by="primerid")
    # colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
    # de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "BH")

    fcHurdle <- merge(summaryDt[contrast=='disease_statusCov19' & component=='H',.(primerid, `Pr(>Chisq)`)],
        summaryDt[contrast=='disease_statusCov19' & component=='logFC', .(primerid, coef, ci.hi, ci.lo, z)], by='primerid')
    
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

    fcHurdleSig <-fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD]

    return(fcHurdleSig)
}


# merge all DE results into a single data.frame 
# this function also cleans the cell type annotations 
merge_organs <- function(tissues = c('liver','lung','kidney','prostate','heart'), 
                            file_id = 'DGE_allcelltypes_ngenes_regression.csv', 
                            output_file_id = '', covariate = 'ngenes'){

        # Read DE output for MAST and merge into a single data.frame
        organs_dge = list()
        for(t in 1:length(tissues)){
            organs_dge[[t]] = read.csv(file =paste0(DGE_DIR, tissues[t], '_', file_id), row.names = 1, header = T)
            organs_dge[[t]]$tissue = tissues[t]
        }
        master_df <- do.call(rbind, organs_dge )

        master_df %>% mutate(full_cell_type = paste0(tissue, ": " , cell_type)) -> master_df
        master_df %>% rename(gene = primerid, pval = fdr, log2fc = coef) -> master_df 
        master_df$method = 'MAST'
        master_df$covariate = covariate 
        master_df %>% dplyr::select(gene, pval, log2fc, cell_type, tissue, method, covariate) -> master_df

        #write.csv(master_df, file = paste0(DGE_DIR, 'all_celltypes_DGE_',output_file_id,'_MAST.csv'), quote = F, row.names = F)
        return(master_df)
}

# Additional DGE methods 

# # # # 
 # # # #
# # # #
# ZINBWAVE 

# zinbwave zero-inflated model 
# this package creates an object that works as wrapper for other DGE packages to work 
# only con is that it takes really long time to run-- make sure to filter the number of genes (3-4k would be fine)
dge_zib_object <- function(core = c() ){
    zinb <- zinbFit(core, X = '~ disease_status', epsilon = 1e12)

    counts = assay(core)
    weights_zinbwave = computeObservationalWeights(zinb, counts)

    # takes a long time to run 
    tissue_zinb <- zinbwave(core, fitted_model = zinb, K = 2, epsilon=1e12,observationalWeights = TRUE)
    return(tissue_zinb)

}

# # DeSeq2 works well with zinbfit (only takes a long time)
# DeSeq2 over a zinb object created by zinbwave. This works fine and returns pvalues and lfcSE (which works with iDEA)
dge_deseq2 <- function(){
    weights <- assay(tissue_zinb, "weights")

    design = model.matrix(~ colData(core)$disease_status)

    dds <- DESeqDataSet(tissue_zinb, design = design)

    dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)

    res <- lfcShrink(dds, coef = 2,
                    type = "normal")

    res %>% as.data.frame  %>% arrange(padj)

    return(res)
}


# # # PREPARE for pathway enrichment analysis 
# convert to iDEA format 
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

# convert p-values from scanpy to beta, beta_var required for iDEA 
scanpy_idea <- function(res = data.frame() ){
    
    pvalue <- res$pvals_adj
    zscore <- qnorm(pvalue/2.0, lower.tail = FALSE)

    beta <- res$logfoldchanges
    se_beta <- abs(beta/zscore)
    beta_var <- se_beta^2
    summary_scanpy  = data.frame(beta = beta, beta_var = beta_var)
    rownames(summary_scanpy) <- res$names
    return(summary_scanpy)

}
# get summary data frame 
#summary_mast <- mast_idea(res= de_res)

# write output 
#write.csv(summary_mast, file =paste0(output_dir, "DGE_", tissue,"_",cell_type,".csv"))


# GENE SET ENRICHMENT 





# GSEA on Reactome pathways 
# this is the main function that takes a single tissue and performs pathway enrichment based on its DGE results
# It will find pathway enrichment per cell type for all cell types in this organ 

# Basic function with Reactome pathways 
gseaReactome <- function(zscores = c()){
	# find pathways for our list of DEG 
	pathways <- reactomePathways(names(zscores))
	# run the pathway enrichment -- really fast 
	fgseaRes <- fgsea(pathways, zscores, maxSize=500, eps = 0 )
	return(list(fgseaRes, pathways)) 
	}


# DE_df is a data.frame of DE results coming from our pipepilines
# standard columns are expected including: pvalue, gene, log2fc, cell_type, tissue 
gsea_tissue <- function(tissue ="liver", 
                        DE_df = data.frame() ){
    # Convert to ENTREZ ID from Symbol 

    library(org.Hs.eg.db)
    library(fgsea)
    library(reactome.db) 

    hs <- org.Hs.eg.db

     # a list whose names are symbols and values are ENTREZ IDs 
    xx <- as.list(org.Hs.egALIAS2EG)
    entrez_db <- xx[!is.na(xx)]


    organ_files = list.files() 
    #aa <- read.csv(input_path, header = T, row.names = 1)
    # read from DGE results -- data.frame with p-values and FC 
    #my.symbols <- row.names(aa)

    organ_gsea = list() 
    # all cell types in this tissue 
    DE_df$cell_type %>% unique -> cell_type_list

    for(c in 1:length(cell_type_list)){
        
        cell_type_dge <- DE_df %>% dplyr::filter(cell_type ==cell_type_list[c])

        my.symbols <- cell_type_dge$gene # as per our data.frame format 


        ## returns all the ENTREZ IDs
        name_map <- mapEntrez(my.symbols, entrez_db) 


        # Filter DGE based on gene name 
        row.names(cell_type_dge) <- cell_type_dge$gene
        dge_fil <- cell_type_dge[name_map$SYMBOL, ] 


        # Compute all zscores
        print( paste('Ranking all genes.. ', cell_type_list[c]) )
        zscores <- make_ranks(dge_fil, name_map) 

        # returns a list 
        print( paste('GSEA.. ', cell_type_list[c]) )
        res_list <- gseaReactome( zscores)
        gsea_res = res_list[[1]]

        export_res <- gsea_res[,1:7] %>% as.data.frame

        gene_strings = c() 
        for(i in 1:dim(export_res)[1]){
            genes<-name_map[which(name_map$ENTREZID %in% gsea_res$leadingEdge[i][[1]]),]$SYMBOL
            gene_strings[i] <- paste(genes, collapse =',')
            
        }
        export_res$genes <- gene_strings

        export_res$tissue = tissue 
        export_res$cell_type = cell_type_list[c]

        organ_gsea[[c]] = export_res 
    }

    export_df <- do.call(rbind, organ_gsea)
    return(export_df)
}

# takes a list of symbols as input
# entrez_map is generated previously with: xx <- as.list(org.Hs.egALIAS2EG)
mapEntrez <- function(my.symbols = c(), entrez_map = c() ){ 
            ENTREZ_ID  = c() 
        for( i in 1:length(my.symbols)){

            entrez_idx <- which(names(xx) == my.symbols[i])
            if(length(entrez_idx)){   
            
                ENTREZ_ID[i] <- xx[[entrez_idx]][1]
            }else{
                ENTREZ_ID[i] = "NA"
            }
        }

        name_map = data.frame(SYMBOL = my.symbols, ENTREZID = ENTREZ_ID)

        name_map %>% dplyr::filter( ENTREZID != "NA") -> name_map
        return(name_map)
}

# compute zscores for every gene basde on their p-value 
make_ranks <- function(dge_fil  = data.frame() , name_map = data.frame() ){ 

        # 1. Converts pval to zscore
        # 2. creates data.frame with log 2 FC and scores, using ENTREZ ID
        pval = dge_fil$pval
        zscores = qnorm(1 - (pval/2))
        zscores[is.infinite(zscores)] = max(zscores[!is.infinite(zscores)])
        logfc = dge_fil$log2fc
        zscores[logfc<0] = -zscores[logfc<0]
        names(zscores) = name_map$ENTREZID
        return(zscores)
}


# merge similar cell types across tissues and find shared signatures of changes in gene expression 
# 1. load all DGE files 
# 2. Create a unique dictionary of cell types 
# 3. GSEA on a cell-type basis across tissues 

gsea_celltype <- function( input_file = paste0(DGE_DIR, 'all_celltypes_DGE_clean_MAST.csv')){



}