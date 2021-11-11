# save counts as .txt 
DATA_DIR = '/mnt/ibm_lg/covid_tissue_atlas/data/tissue_objects/all_tissues/'
OUT_DIR = '/mnt/ibm_lg/covid_tissue_atlas/notebooks/cellphone/'

setwd(OUT_DIR) 

tissue = 'lung'



load(paste0(DATA_DIR, tissue, "/r_files/seurat_obj.rds") ) 

# split object into disease and healthy 
obj_list = SplitObject(seurat_obj, split.by ='disease_status')

sample_names = names(obj_list) # condition 
	# let's make new colnames to avoid problems 

for(p in 1:length(sample_names) ){

	colnames_ = colnames(obj_list[[p]][['RNA']]@data)
	colnames_ <- paste('X', colnames_, sep = "")
	colnames_ = c('Gene',  colnames_)
	

# Prepare meta data files 
	genes = row.names(obj_list[[p]][['RNA']]@data)

	write.table(data.frame('Gene' = genes, obj_list[[p]][['RNA']]@data), 
		quote = F, row.names = F, sep = '\t', col.names = colnames_, 
		file = paste0(OUT_DIR, tissue,'/counts_',sample_names[p] , '.txt'))

	
	# save meta file 
	# use the modified cell ids 
	meta<- data.frame('Cell' = colnames_[-1] , 
										'cell_type' = paste( obj_list[[p]]@meta.data %>% pull(cell_type_annotation), obj_list[[p]]@meta.data %>% pull(disease_status) ) )
	
	write.table(meta, , quote = F, row.names = F, 
		sep = '\t', file = paste0(OUT_DIR, tissue, '/meta_',sample_names[p] , '.txt'))

}

setwd( paste0(OUT_DIR,tissue) ) 
for(p in 1:length(sample_names) ){
#Run CELLPHONE 
	# to run cellphone DB
  
	# run CellPhone DB
	phone_cmd <- paste0( 'cellphonedb method statistical_analysis ',
											  'meta_', sample_names[p], '.txt counts_', sample_names[p], '.txt --counts-data gene_name ',
										 '--output-path out_',sample_names[p]  )
  system(phone_cmd) 
}
setwd(OUT_DIR) 

# works only for one condition 
makeCellPairs <- function(cell_types= c('macrophage', 'endothelial'), 
													sample_name = "Control" , 
													sep1 = ".", 
													sep2 = "."){
	# concat with .
	x <- paste(cell_types, sample_name, sep= sep1)
	# all pairs
	all_pairs <- permutations(n=2 ,r = 2,  v = x, repeats.allowed = T)
	# final list of pairs
	Sig_columnpair <- apply(all_pairs, 1, function(x){paste(x, collapse = sep2)} )
  return(Sig_columnpair)
}
	

library(gtools)

this_condition = "Control"
setwd(paste0(OUT_DIR, tissue ) ) 

Sig_columnpair = makeCellPairs(cell_types = c('macrophage', 'endothelial'), 
														   sep1 =".", sep2=".", 
															 sample_name = this_condition) 

#EXPORT FILES 
out_path = paste0("out_",this_condition, "/")

DF<-read.table(paste0(out_path, "pvalues.txt"),sep="\t",header = T)

rownames(DF)<-make.unique(as.vector(DF$interacting_pair))
DF_orig<-DF

DF<-DF[,colnames(DF) %in% Sig_columnpair]
Sig_rows<-NULL

for (i in 1:length(Sig_columnpair)) {
  Sig_rows_temp<-rownames(DF)[DF[,i]<0.05]
  if (length(Sig_rows_temp)>=1) {
  Sig_rows<-union(Sig_rows,Sig_rows_temp)}
}

DF_orig<-DF_orig[rownames(DF_orig) %in% Sig_rows,]

#write all files 
Sig_columnpair = makeCellPairs(cell_types = c('macrophage', 'endothelial'), 
														   sep1 =" ", sep2="|", 
															 sample_name = this_condition) 

write.table(Sig_columnpair,paste0(out_path, "columns_direction1.txt"),sep="\n",col.names = F,row.names = F,quote = F)


write.table(DF_orig,paste0(out_path, "Filtered_pvalues_direction1.txt"),sep="\t",row.names = T,col.names = T,quote = F)
write.table(Sig_rows,paste0(out_path,"rows_direction1.txt") ,sep="\n",col.names = F,row.names = F,quote = F)



# NOW we can run  CPDB

# cellphonedb plot dot_plot --rows out_Cov19/rows_direction1.txt  --columns out_Cov19/columns_direction1.txt --output-path=out_Cov19 --pvalues-path=out_Cov19/pvalues.txt --means-path=out_Cov19/means.txt
# cellphonedb plot dot_plot --rows out_Control/rows_direction1.txt  --columns out_Control/columns_direction1.txt --output-path=out_Control --pvalues-path=out_Control/pvalues.txt --means-path=out_Control/means.txt