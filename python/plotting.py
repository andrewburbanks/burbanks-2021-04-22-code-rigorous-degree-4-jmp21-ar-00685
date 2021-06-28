import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')

nice_fonts = {
    # Use LaTex to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 0.75,
    "savefig.format": 'pdf',
    "figure.dpi": 180, # larger values make the image larger in the jupyter notebook
    "savefig.dpi": 300,
    "path.simplify": True,
    'legend.numpoints': 1,
    'legend.frameon': False,
    'legend.handletextpad': 0.5,
}

def title(txt):
    print(f'Figure: {txt}')

def init_plotting(plt):
    plt.rcParams.update(nice_fonts)
    plt.title = title

def set_size(width, fraction=1, subplots=[1,1], ratio=None):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Parameters
    ----------
    width: float
            Width in pts
    fraction: float
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    
    cms_per_inch = 2.54
    inches_per_pt = 1 / 72.27
    
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    elif width == 'pnas':
        width_pt = 246.09686
    elif width == 'prl':
        # single column of two-column prl
        # 244.69 #print(width_pt)
        width_pt = 8.6/cms_per_inch/inches_per_pt
    else:
        width_pt = width
        
    # Width of figure
    fig_width_pt = width_pt * fraction

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2

    if ratio is None:
        ratio = golden_ratio
    
    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * ratio * (subplots[0] / subplots[1])

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim

def savefig(plt, filename, **kwargs):
    plt.savefig(filename, bbox_inches='tight', dpi=300, **kwargs) #, transparent=True)

