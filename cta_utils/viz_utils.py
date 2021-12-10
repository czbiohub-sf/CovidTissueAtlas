import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys



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

    sc.settings.verbosity = 3            
    sc.set_figure_params(dpi=150)
    sc.settings.figdir = FIG_DIR


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


def stylize_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top='off', direction='out', width=0)
    ax.yaxis.set_tick_params(right='off', direction='out', width=0)
    
def custom_barchart(ax, x, y, error, xlims, ylims, error_kw, color='lightblue', width=0.75):
    """Customized bar chart with positive error bars only."""

    error = [np.zeros(len(error)), error]

    ax.bar(x, y, color=color, width=width, yerr=error, error_kw=error_kw, align='center')

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    
    
def stylize_axes(ax, title, xlabel, ylabel ): #, xticks, yticks, xticklabels, yticklabels):
    """Customize axes spines, title, labels, ticks, and ticklabels."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top='off', direction='out', width=0)
    ax.yaxis.set_tick_params(right='off', direction='out', width=0)
    
    ax.set_title(title)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
#     ax.set_xticks(xticks)
#     ax.set_yticks(yticks)
    
#     ax.set_xticklabels(xticklabels)
#     ax.set_yticklabels(yticklabels)
# Vizualiation

def cell_type_barplot(adata_obj, cell_type_label='annotations_v1',group_by='disease_status', cluster_lab  ='leiden'):
    cell_type_counts = adata_obj.obs.groupby(by=[cell_type_label,group_by]).count()[cluster_lab].reset_index()

    ax = sns.histplot(
        cell_type_counts,
        y=cell_type_label,
        # Use the value variable here to turn histogram counts into weighted
        # values.
        weights=cluster_lab,
        hue=group_by,
        multiple='dodge',
        palette=['#24b1d1', '#ae24d1'],
        # Add white borders to the bars.
        edgecolor='white',
        # Shrink the bars a bit so they don't touch.
        shrink=0.8
    )

