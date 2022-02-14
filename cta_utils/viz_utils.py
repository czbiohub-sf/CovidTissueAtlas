import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import scanpy as sc



# we use this header in all notebooks 

def setup_fig_params(this_sample =""): 
    # directory final object and figures will be saved to
    FIG_DIR = '/mnt/ibm_lg/covid_tissue_atlas/figures/figure2/' + this_sample + '/'

    # if directory doesn't exist, make directory
    if not os.path.exists(FIG_DIR):
        os.makedirs(FIG_DIR)

    #print package versions
    print('package versions:')
    print('\n'.join(f'{m.__name__} {m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))

    # make sure you use these two parameters before saving figures to pdf
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Arial'

    return FIG_DIR 



# Plotting style function (run this before plotting the final figure)
def set_plotting_style():
    plt.style.use('seaborn-paper')
    plt.rc('axes', labelsize=12)
    plt.rc('axes', titlesize=12)
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('legend', fontsize=10)
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('xtick.major', pad=2)
    plt.rc('ytick.major', pad=2)
    plt.rc('mathtext', fontset='stixsans', sf='sansserif')
    plt.rc('figure', figsize=[10,9])
    plt.rc('svg', fonttype='none')


# def stylize_axes(ax):
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)

#     ax.xaxis.set_tick_params(top='off', direction='out', width=0)
#     ax.yaxis.set_tick_params(right='off', direction='out', width=0)
    
def custom_barchart(ax, x, y, error, xlims, ylims, error_kw, color='lightblue', width=0.75):
    """Customized bar chart with positive error bars only."""

    error = [np.zeros(len(error)), error]

    ax.bar(x, y, color=color, width=width, yerr=error, error_kw=error_kw, align='center')

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    
# we use this in the main notebooks
def stylize_axes(ax, title, xlabel, ylabel ): #, xticks, yticks, xticklabels, yticklabels):
    """Customize axes spines, title, labels, ticks, and ticklabels."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top='off', direction='out', width=0)
    ax.yaxis.set_tick_params(right='off', direction='out', width=0)
    
    ax.set_title(title)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    

# Vizualiation

def cellType_barplot(adata = [], ax1 = [], 
                     use_annotation = 'short_cell_type',
                     condition = 'disease_status'): 

    stacked_df = adata.obs[[condition, use_annotation]]
    cluster_disease = stacked_df.pivot_table(index=use_annotation, columns=[condition], aggfunc='size')
    cluster_disease_pct = cluster_disease.div(cluster_disease.sum(axis=0), axis=1) * 100

    n_types = len(set(adata.obs.cell_type_annotation))
    cluster_disease_pct.plot(kind='bar', stacked=False, 
                            grid = False, cmap = 'Spectral' , 
                            ax = ax1)

    ax1.legend().set_visible(False) # we don't need the legend for the integrated figure    
    ax1.set_ylabel('Proportion (%)');
    ax1.set_xlabel('Cell type');
    return ax1


def pretty_umap(adata, group_by = 'condition', use_pal='Spectral',


                use_title= "COVID status", 
                ax1 =[]): 

    sc.pl.umap(adata, color=group_by, add_outline=True, legend_loc='right margin',
                   legend_fontsize=10, legend_fontoutline=2,frameon=False,
                   title= use_title, palette=use_pal,
                   ax=ax1, show = False)

    return ax1