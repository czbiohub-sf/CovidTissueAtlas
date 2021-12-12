import anndata
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import os
import scvi

# Data location -- final h5ads 
DATA_DIR = '/mnt/ibm_lg/covid_tissue_atlas/data/tissue_objects/all_tissues/'
# File names for each tissue
# keep this dictionary updated!
sample_dict = {'liver':  'liver_CTA_only.h5ad',
               'lung':   'lung_CTA_TS_integrated.h5ad',
               'kidney': 'kidney_CTA_KA_integrated.h5ad',
               'heart': 'heart_CTA_only.h5ad',
               'prostate': 'prostate_CTA_only.h5ad'}


# COVID genes of interest
# covid entry factors 
def get_covid_genes(): 
    return {'Entry genes' : ['ACE2', 'TMPRSS2', 'CTSL'],
    'Entry receptors' : ['ACE2', 'BSG', 'NRP1', 'DPP4', 'ANPEP', 'ENPEP', 'AGTR2'], 
    'Entry proteases' : ['TMPRSS2', 'CTSL', 'FURIN', 'PRSS1'] } 

# Function loading the processed datasets for figure making 
def load_final_tissue(this_sample, cta_only = True): 
    # Load final object (subset from the CTA atlas)
    if(cta_only): 
        sample_name  = this_sample + "_CTA_only.h5ad"
    else: 
        sample_name = sample_dict[this_sample] # Load an integrated object containing external data (for DE )

    adata_path = DATA_DIR + this_sample + '/'  + sample_name 
    
    print("Loading " + this_sample +  " dataset located at " + adata_path ) 
    

    adata = sc.read_h5ad(adata_path)

    # convet cell_ontology_class to short name 
    short_names =short_celltype_names(this_sample)
    short_cell_names = [short_names[e] for e in adata.obs.cell_ontology_class  ]
    adata.obs['short_cell_type'] = short_cell_names

    return adata 




def read_cellbender_h5ad_to_combined_adata(data_path):
    '''
    Reads in a data path to directory containing cellbender-processed files to be combined into adata object. Returns one adata object
    containing all h5ad objects in data path. No filtering is performed on adata object.
    '''
    adata = sc.AnnData()

    data_folders = [f for f in os.listdir(data_path) if not f.startswith('.')]

    for d in data_folders:
        print(d)
        adataaux = sc.read_10x_h5(data_path + d)
        bfcmingenes = adataaux.shape[0]

        print(d+ " size")
        display(adataaux.shape)

        adataaux.obs['10X_run'] = d.split('_')[1]


        try:
            adataaux.obs['disease_status'] = d.split('_')[0]
            adataaux.obs['lab'] = d.split('_')[1]
            adataaux.obs['tissue'] = (str.lower(d.split('_')[2])).split('.')[0]
            adataaux.obs['Barcode'] = adataaux.obs_names 
            adataaux.obs_names = adataaux.obs['10X_run'] + '_' + adataaux.obs_names

        except:
            adataaux.obs['disease_status'] = d.split('-')[0]
            adataaux.obs['lab'] = d.split('-')[1]
            adataaux.obs['tissue'] = str.lower(d.split('-')[2].split('_')[0])


        try:
            adataaux.obs['replicate'] = d.split('_')[4]
        except:
            adataaux.obs['replicate'] = 1

        try:
            adataaux.var_names_make_unique()
            adata = adata.concatenate(adataaux,join='inner', index_unique = None)
        except:
            adata = adataaux.copy()

    adata.raw = adata.copy()

    adata_gencode = adata.copy()

    adata.var['gene_symbol'] = [g.split("_")[-1] for g in adata.var.index]
    adata.var = adata.var.set_index('gene_symbol')

    return adata

def read_decontx_h5ad_to_combined_adata(data_path):
    '''
    Reads in a data path to directory containing 10x files to be combined into adata object. Returns one adata object
    containing all 10x files in data path. No filtering is performed on adata object.
    '''
    adata = sc.AnnData()

    data_folders = [f for f in os.listdir(data_path) if not f.startswith('.')]

    for d in data_folders:
        print(d)
        #runs = os.listdir(data_path+d)
        #runs = [f for f in os.listdir(data_path + d) if not f.startswith('.')]

        #print(r)
        adataaux = sc.read_h5ad(data_path + d)
        bfcmingenes = adataaux.shape[0]

        print(d+ " size")
        display(adataaux.shape)

        adataaux.obs['10X_run'] = d.split('.')[0]


        try:
            #adataaux.obs['sample_type'] = d.split('_')[3]
            adataaux.obs['disease_status'] = d.split('_')[0]
            adataaux.obs['lab'] = d.split('_')[1]
            adataaux.obs['tissue'] = (str.lower(d.split('_')[2])).split('.')[0]
            adataaux.obs_names = adataaux.obs['10X_run'] + '-' + adataaux.obs['Barcode']

        except:
            #adataaux.obs['sample_type'] = d.split('-')[2].split('_')[1]
            adataaux.obs['disease_status'] = d.split('-')[0]
            adataaux.obs['lab'] = d.split('-')[1]
            adataaux.obs['tissue'] = str.lower(d.split('-')[2].split('_')[0])


        try:
            adataaux.obs['replicate'] = d.split('_')[4]
        except:
            adataaux.obs['replicate'] = 1

        try:
            adata = adata.concatenate(adataaux,join='outer', index_unique = None)
            #adata.obs = adata.obs.drop('batch',axis=1)
        except:
            adata = adataaux.copy()

    adata.raw = adata.copy()

    #set adata.X as decontX counts
    adata.X = adata.layers['decontXcounts']

    adata_gencode = adata.copy()

    adata.var['Symbol'] = [g.split("_")[1] for g in adata.var['Symbol']]
    adata.var = adata.var.set_index('Symbol')

    adata.obs['decontX_contamination'] = adata.obs['decontX_contamination'].astype(str)
    adata.obs['decontX_clusters'] = adata.obs['decontX_clusters'].astype(str)
    return adata

def filter_adata(adata, min_n_counts, min_n_genes, filter, mt_cutoff = 50, filter_max = 0):

    adata.var_names_make_unique()
    adata_no_mt = adata[:, ~adata.var_names.str.startswith('MT-')].copy()

    # filter genes to speed up SCVI
    sc.pp.filter_genes(adata_no_mt, min_cells=3)
    
    # filter cells using the no-mt data
    sc.pp.filter_cells(adata_no_mt, min_counts= min_n_counts)

    if(filter_max):
        sc.pp.filter_cells(adata_no_mt,max_counts = filter_max)

    sc.pp.filter_cells(adata_no_mt, min_genes = min_n_genes )

    adata = adata[adata.obs.index.isin(adata_no_mt.obs.index)].copy()

    if filter == 'normal':
        return adata

    elif filter == 'percent_mito':
        #remove cells with greater than mt_cutoff % mito genes
        adata = adata[adata.obs.pct_counts_mt < mt_cutoff, :]
        return adata

    elif filter == 'remove_all_mito_genes':
        return adata_no_mt
    
def run_scvi(adata, n_layers, n_latent, batch_label ='10X_run'):
    """
    Takes in adata with raw counts. Saves raw counts in a separate layer. Perform normalization with standard pipeline,
    and prepares SCVI model.Train SCVI model. Input number of layers (n_layers) and latent spaces
    (n_latent) to use.
    """
    #SCVI preprocessing
    # normalize data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep full dimension safe

    # For batch correction: Filter samples that have less than 200 cells
    sample_id = batch_label # decontX object
    sample_size = adata.obs[sample_id].value_counts()
    large_samples = sample_size[sample_size>200].reset_index()['index'].astype(str)

    adata = adata[adata.obs[sample_id].isin(large_samples)].copy()

    # Train SCVI
    scvi.data.setup_anndata(adata, layer="counts", batch_key=batch_label)
    vae = scvi.model.SCVI(adata, n_layers=n_layers, n_latent=n_latent) # or even 50

    # takes about 40 min
    vae.train()

    adata.obsm["X_scVI"] = vae.get_latent_representation() # latent vectors scvi, n = 30
    return adata


def compute_neighbors(adata, use_rep, k = 15 , use_metric = 'cosine'):
    """
    Compute neighbors using decontx (decontX_UMAP) or SCVI (X_scVI) -corrected obsm.
    """
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors = k , metric = use_metric )
    sc.tl.leiden(adata)
    #sc.tl.louvain(adata)
    sc.tl.umap(adata)
    return adata

def read_h5ad_to_combined_BroadData(data_path, mode ='h5ad'):
    '''
    Reads in a data path to directory containing 10x files to be combined into adata object. Returns one adata object
    containing all 10x files in data path. No filtering is performed on adata object.
    '''
    adata = sc.AnnData()

    data_folders = [f for f in os.listdir(data_path) if not f.startswith('.')]

    for d in data_folders:
        print(d)

        if(mode=='h5ad'):
            adataaux = sc.read_h5ad(data_path + d)
        else:
            adataaux = sc.read_10x_h5(data_path + d)

        # set gene index as gene Symbol
        adataaux.var['Symbol'] = adataaux.var.index
        adataaux.obs['Barcode'] = adataaux.obs.index

        bfcmingenes = adataaux.shape[0]

        print(d+ " size")
        display(adataaux.shape)

        # Donor ID
        adataaux.obs['10X_run'] = d.split('-')[2] #d.split('.')[0]
        adataaux.obs['lab'] = d.split('_')[0] #unique sample ID
        adataaux.obs_names= [ d.split('_')[1] + '-' + b.split('-')[0]   for b in adataaux.obs['Barcode']] #adataaux.obs['cell_id']
        adataaux.var_names_make_unique()

        try:
            adata = adata.concatenate(adataaux,join='outer', index_unique = None)
            #adata.obs = adata.obs.drop('batch',axis=1)
        except:
            adata = adataaux.copy()


    adata.raw = adata.copy()
    adata_gencode = adata.copy()

    adata.var['Symbol'] = [g.split("_")[1] for g in adata.var['Symbol']]
    adata.var = adata.var.set_index('Symbol')


    return adata

def concat_2anndata(adata1,adata2):
    min_cols = list(set(adata1.obs.columns) & set(adata2.obs.columns))
    min_vars = list(set(adata1.var.columns) & set(adata2.var.columns))

    if(adata1.layers):
        adata = anndata.AnnData(X = adata1.layers['counts'],
                       obs= adata1.obs[min_cols],
                       var= adata1.var[min_vars])
    else:
        adata = anndata.AnnData(X = adata1.X,
                       obs= adata1.obs[min_cols],
                       var= adata1.var[min_vars])

    adataaux = anndata.AnnData(X = adata2.layers['counts'],
                              obs = adata2.obs[min_cols],
                              var = adata2.var[min_vars])

    adata = adata.concatenate(adataaux, join='inner',index_unique = None)

    adata_com = adata1.concatenate(adata2, join='inner',index_unique = None)

    return adata


def get_markers(tissue =""):
    all_markers = {
        'prostate': {'luminal':         ['KLK3','NKX3-1','ACPP','KRT8'],
                     'basal':           ['KRT5','NOTCH4','TP63','DST'],
                     'endothelial':     ['VWF','PECAM1','SELE'],
                     'smooth muscle':   ['ACTA2', 'MYH11', 'TPM2'],
                     'fibroblast':      ['DCN', 'C7', 'PDGFRA']},
        'liver':    {'hepatocyte':      ['ALB', 'APOA1', 'HAL', 'SDS','HAMP'],
                    'stromal':          ['ACTA2', 'COL1A1', 'DCN', 'COLEC11', 'SPARC'],
                    'endothelial':      ['PECAM1', 'CLEC4G', 'CLEC4M', 'FLT4','FCGR2B'],
                    'smooth muscle':    ['ACTA2', 'MYH11', 'CNN1', 'TAGLN'],
                    'kupffer':          ['CLEC4F', 'CD5L', 'VSIG4',  'MST1R', 'MAFB'],
                    'NK cell':          ['PTPRC', 'KLRB1', 'KLRF1', 'GZMA']},
        'heart':    {'endothelial':     ['VWF', 'PECAM1', 'CDH5'],
                    'cardiac muscle':   ['TNNT2', 'TNNI1', 'TNNI3', 'ACTC1'],
                    'fibroblast':       ['DCN', 'LUM', 'PDGFRA'],
                    'smooth muscle':    ['MYH11', 'TAGLN', 'ACTA2', 'CNN1', 'RGS5', 'AGT', 'CSPG4'],
                    'macrophage':       ['CD163', 'PTPRC', 'MARCO', 'FCGR3A'],
                    'T-cell CD8':       ['CD8A', 'IL7R', 'PTPRC']},
        'lung': {   'neuron':           ['CACNA1A', 'NRXN1', 'NRXN3', 'CNTNAP2', 'RBFOX1'],
                    'AT-I':                ['AGER', 'PDPN', 'CLIC5'],
                    'AT-II':               [ 'SFTPB', 'SFTPC', 'SFTPD', 'MUC1', 'ETV5'],
                    'basal cell':          ['KRT5', 'KRT14', 'TP63', 'DAPL1'],
                    'artery cell':         ['GJA5'],
                    'macrophage':          ['MARCO', 'MSR1', 'MRC1'],
                    'fibroblast':          ['COL1A1', 'PDGFRA'],
                    'NK cell':             ['KLRD1', 'NKG7', 'TYROBP'],
                    'macrophage':          ['MARCO', 'MSR1', 'MRC1']},
        'kidney':{  'p conv tubule':      ['CUBN', 'SLC4A4', 'SLC34A1', 'MIOX'],
                    'alpha':              [ 'SLC4A1','OXGR1','ATP6AP2'],
                    'asc loop Henle':     ['SPP1', 'SLC12A1', 'UMOD'],
                    'd conv tubule':      ['SLC12A3', 'WNK1', 'WNK4', 'SLC12A1'] ,
                    'B cell' :            ['CD79A', 'CD24', 'MS4A1', 'CD19'] ,
                    'podocyte':           ['TJP1', 'WT1', 'NPHS1', 'NPHS2', 'THSD7A', 'SYNPO', 'CD2AP', 'COL4A3', 'GOLIM4', 'PODXL', 'PTPRO', 'FOXD1'],
                    'Endothelial':        ['GJA4','BMX','ACKR1','AQP1','PROX1', 'FLT4', 'CA4','PDPN','PECAM1','CD34','VWF']}
    }
    return all_markers[tissue]

# for a given organ it returns a dictionary with short celltype names 
# it should be then used as pd.rename(columns = dict)
# or in foor loop to create a short_name column in adata.obs 
def short_celltype_names(tissue =""): 
    short_names = {
            'lung': {
                'type i pneumocyte': 'AT-I', 
                'type ii pneumocyte' : 'AT-II', 
                'lung ciliated cell': 'Ciliated',
                't cell' : 'T-cell',
                'lung macrophage': 'Macrophage',
                'lung endothelial cell': 'Endothelial',
                'fibroblast of lung': 'Fibroblast',
                'dendritic cell': 'Dendritic',
                'basal cell': 'Basal',
                'b cell': 'B cell'





            }, 
            'kidney': {}
    
    }
    return short_names[tissue]

