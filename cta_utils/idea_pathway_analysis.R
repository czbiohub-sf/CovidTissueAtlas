# pathway (gene set) enrichemnt analysis
# input for iDEA: data frame col1=log_fold_change , col2=p.value
# depending on the DGE method used, we create a table

library(doParallel)
library(BiocParallel)
library(iDEA)

input_file = ''

output_dir = '/mnt/ibm_lg/covid_tissue_atlas/results/DGE/'
input_file = 'DGE_liver_hepatocyte.csv'


NCORES = 25 
registerDoParallel(NCORES)
register(DoparParam())



summary_data = read.csv( paste0(output_dir,input_file) , header = T, row.names = 1)

data('humanGeneSets')
annotation_data = humanGeneSets




genes_common = row.names(summary_data)[which(row.names(summary_data) %in% row.names(annotation_data))]


# make model and fit 
idea <- CreateiDEAObject(summary_data[genes_common,], 
                        annotation_data[genes_common,], 
                        max_var_beta = 100, 
                        min_precent_annot = 0.0025, num_core=50)

idea <- iDEA.fit(idea,
                 fit_noGS=FALSE,
                 init_beta=NULL, 
                 init_tau=c(-2,0.5),
                 min_degene=5,
                 em_iter=15,
                 mcmc_iter=1000, 
                 fit.tol=1e-5,
                 modelVariant = F,
                 verbose=TRUE)
# correct p-values 
idea<-iDEA.louis(idea)