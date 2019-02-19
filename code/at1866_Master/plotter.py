"""plotting functionality
"""
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

try:
    from IPython.display import set_matplotlib_formats
    set_matplotlib_formats('retina')
except ImportError:
    pass



def plotter(fig=None, base=None, plot_aux=None, buffered=None, model=None,
            pt1_size=None, pt2_size=None, plot_res=None, save_fig=False,
            title=None, area=None, census_geo=None, sa_style='spider_lines',
            figsize=(10,10)):
    """ Top-level scenario plotter for location analytics.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        complete figure to plot. Default is None.
    base : matplotlib.axes._subplots.AxesSubplot
        individual axis to plot. Default is None.
    plot_aux : dict
        model data parameters dataframes to plot keyed by
        descriptive names. Default is None.
    plot_res : dict
        model data results dataframes to plot keyed by
        descriptive names. Default is None.
    buffered : see
        buffer distance from roads segments in `plot_base`. Default is None.
    census_geo : str
        spatial unit of census geography.
    pt1_size : float or float
        size of points to plot. `pt1_size` should always be the
        larger between `pt2_size` and `pt1_size`. Default is None.
    pt2_size : float or float
        size of points to plot. Default is None.
    model : pulp.pulp.LpProblem
        location model
    area : str
        location of model
    title : str
        title. Default is None.
    sa_styleplot : str
        method for plotting service area. Default is spider_lines.
        Option is concave_hull, or dict for network area.
    figsize : tuple
        Figure size for plot. Default is (10,10).
    save_fig : bool
        Default is False.
    
    Returns
    -------
    add_to_legend : list
        items to add to legend
    """
    if area == 'Leon_FL':
        figsize = (25,20)
    
    for_multiplot = True
    if not fig and not base:
        for_multiplot = False
        fig, base = plt.subplots(1, 1, figsize=figsize)
    
    # add title
    if not for_multiplot:
        if model:
            title += ' - ' + model.name
        base.set_title(title, size=20)
    else:
        base.set_title(model.name, size=20)
    
    # plot non-results data
    if plot_aux:
        for k, df in plot_aux.items():
            if k == 'census_polys':
                df.plot(ax=base, color='b', alpha=.035, linestyle='--',
                        linewidth=2, zorder=1)
            if k == 'streets':
                street_width = 2
                if type(sa_style) == dict:
                    street_width = .25
                df.plot(ax=base, lw=street_width, color='k', zorder=1)
            if k == 'buffer':
                df.plot(ax=base, color='y', lw=.25, alpha=.25, zorder=1)
            if k == 'cli_tru':
                if plot_res:
                    df = df[df[model.name+'_sol'] == 'closed']
                    psize = pt2_size/6.
                    pcolor = 'k'
                else:   
                    n_cli = df.shape[0]
                    psize = pt1_size
                    pcolor = 'b'
                df.plot(ax=base, markersize=psize,
                        edgecolor='k', color=pcolor)
            if k == 'fac_tru':
                if plot_res:
                    df = df[df[model.name+'_sol'] == 'closed']
                    psize = pt2_size
                    pcolor = 'k'
                    pmarker = '*'
                else:
                    n_cli = df.shape[0]
                    psize = pt1_size
                    pcolor = 'r'
                    pmarker = 'o'
                df.plot(ax=base, markersize=psize,
                        edgecolor='k', color=pcolor,
                        marker=pmarker)
                n_fac = df.shape[0]
            if k == 'cli_snp':
                df.plot(ax=base, markersize=pt2_size,
                        edgecolor='k', color='b', alpha=.75)
            if k == 'fac_snp':
                df.plot(ax=base, markersize=pt2_size,
                        edgecolor='k', color='r', alpha=.75)
        add_to_legend = list(plot_aux.keys())
    else:
        add_to_legend = None
    
    # plot results data
    if plot_res:
        dv_colors = dv_colorset(plot_res['fac_var'].desc_var)
        # facilities
        df = plot_res['fac_var'][plot_res['fac_var']\
                                [model.name+'_sol'] != 'closed']
        alpha = .5
        # decision variable info for legend
        dvs_to_leg = {}
        # plot facilities
        for desc_var in df.desc_var:
            fac = df[df.desc_var == desc_var]
            fac.plot(ax=base, marker='*', markersize=pt1_size*3.,
                     alpha=.8, zorder=3, edgecolor='k',
                     color=dv_colors[desc_var])
            # update decision variable info with set color
            dvs_to_leg[desc_var] = {'color':dv_colors[desc_var]}
        # plot clients & service areas
        for f, c in model.fac2cli.items():
            fc = plot_res['cli_var'][plot_res['cli_var'].desc_var.isin(c)]
            fc.plot(ax=base, markersize=50, edgecolor='k',
                    color=dv_colors[f], alpha=alpha, zorder=2)
            # update decision variable info with set client counts
            dvs_to_leg[f].update({'clients': fc.shape[0]})
            # create service area object
            service_area = service_area_object(df, fc, f, style=sa_style)
            if type(sa_style) == str:
                service_area.plot(ax=base, edgecolor='k', alpha=.2,
                                  color=dv_colors[f], zorder=1)
            else:
                service_area.plot(ax=base, alpha=.2, zorder=1,
                                  color=dv_colors[f], linewidth=10)
    
    else:
        dvs_to_leg = None
    
    # create a shell class to represent FacilityLocationModel if not present
    if not model:
        try:
            model = _ShellModel(plot_aux,
                                ['cli_tru', 'fac_tru', 'census_polys'])
        except (TypeError, KeyError):
            model = None
    # if FacilityLocationModel present add extra attributes for patch creation
    else:
        try:
            key = 'census_polys'
            attr_name = 'n_' + key[:3]
            setattr(model, attr_name, plot_aux[key].shape[0])
        except KeyError:
            pass
    
    if not for_multiplot:
        # create legend patches
        patches = create_patches(model=model, for_multiplot=for_multiplot,
                                 pt1_size=pt1_size, pt2_size=pt2_size,
                                 buffered=buffered, legend_aux=add_to_legend,
                                 dvs_to_leg=dvs_to_leg, census_geo=census_geo)
        add_legend(patches, for_multiplot=for_multiplot)
    add_north_arrow(base, area=area)
    add_scale(base, area=area, for_multiplot=for_multiplot)
    
    if save_fig:
        plt.savefig(model.name+'.png')
    
    # if for a multiplot explicityly return items to add to legend
    if for_multiplot:
        return add_to_legend


class _ShellModel:
    """object to mimic `model` when not present
    """
    def __init__(self, plot_aux, keys):
        for key in keys:
            attr_name = 'n_' + key[:3]
            try:
                setattr(self, attr_name, plot_aux[key].shape[0])
            except KeyError:
                pass


def multi_plotter(models, plot_aux=None, plot_res=None, select=None,
                  title=None, area=None, census_geo=None, net=None, 
                  sa_style='spider_lines', figsize=(14,14), shape=(2,2)):
    """plot multiple base axes as one figure
    
    Parameters
    ----------
    models : list
        solved model objects
    select : dict
        facility-to-selection count lookup.
    shape : tuple
        dimension for subplot array. Default is (2,2).s
    sa_style : str
        see plotter(). Default is 'spider_lines'.
    net : SpaghettiNetwork
        used for create of network service area
    census_geo : see plotter()
    plot_aux : see plotter()
    plot_res : see plotter()
    title : see plotter()
    area : see plotter()
    figsize : see plotter()
    """
    pt1_size, pt2_size = 300, 60
    # convert list of models to array 
    mdls = np.array(models).reshape(shape)
    fig, axarr = plt.subplots(mdls.shape[0], mdls.shape[1],
                              figsize=figsize,
                                      sharex='col',
                                      sharey='row')
    # add super title to subplot array
    plt.suptitle(title, fontsize=30)
    fig.subplots_adjust(hspace=0.1, wspace=0.005, top=.925)
    # create each subplot
    for i in range(mdls.shape[0]):
        for j in range(mdls.shape[1]):
            if net:
                sa_style = {'model':mdls[i,j], 'net':net}
            add_to_legend = plotter(base=axarr[i,j],
                                    plot_aux=plot_aux,
                                    plot_res=plot_res,
                                    model=mdls[i,j],
                                    pt1_size=pt1_size,
                                    pt2_size=pt2_size,
                                    area=area, sa_style=sa_style)
            axarr[i,j].set_aspect('equal')
    add_to_legend = set(add_to_legend)
    # decision variable color set
    dv_colors = dv_colorset(plot_res['fac_var'].desc_var)
    dvs_to_leg = {f: dv_colors[f] for m in models\
                                  for f in m.fac2cli.keys()}
    # set ordered dict of {iloc:fac_var, color, x-selected}
    # *** models[0] can be any of the solved models
    dvs_to_leg = {models[0].fac2iloc[k]:(k,v, select[k])\
                  for k, v in dvs_to_leg.items()}
    dvs_to_leg = OrderedDict(sorted(dvs_to_leg.items()))
    
    try:
        _shell_model = _ShellModel(plot_aux, ['census_polys'])
    except (TypeError, KeyError):
        _shell_model = None
    
    # create legend patches
    patches = create_patches(model=_shell_model,
                             pt1_size=pt1_size, pt2_size=pt2_size,
                             legend_aux=add_to_legend,
                             dvs_to_leg=dvs_to_leg,
                             census_geo=census_geo,
                             for_multiplot=True)
    add_legend(patches, for_multiplot=True)


def add_north_arrow(base, area=None):
    """add a north arrow to an axes
    Parameters
    ----------
    base : see plotter()
    """
    if area == 'Phoenix_grid':
        x, y = 221200, 267200
    elif area == 'Leon_FL':
        x, y = 580000,180000
    elif area == 'Test_Waverly_Leon_FL':
        x, y = 621250, 166500
    elif area == 'Test_Grid_Leon_FL' or area == 'Grid':
        x, y = 619475, 158225
    elif area == 'Test_Small_Leon_FL':
        x, y = 621100, 156125

    arw = 'rarrow, pad=0.25'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    base.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)


def add_scale(base, area=None, for_multiplot=False):
    """add a scale arrow to an axes
    Parameters
    ----------
    base : see plotter()
    area : see plotter()
    for_multiplot : see plotter()
    """
    # set scale anchor
    if area in ['Test_Small_Leon_FL', 'Test_Grid_Leon_FL', 'Small', 'Grid']:
        offset = 15
    else:
        offset = 75
    x, y = base.get_xlim()[0]+offset, base.get_ylim()[0]+offset
    
    # set scale distance and units
    if area in ['Test_Small_Leon_FL', 'Test_Grid_Leon_FL', 'Small', 'Grid']:
        distance, units = 50, 'm'
    
    elif area == 'Test_Waverly_Leon_FL':
        distance, units = .5, 'km'
    
    elif area == 'Leon_FL':
        distance, units = 10, 'km'
        
    elif area == 'Phoenix_grid':
        distance, units = .25, 'km'
        
    '''
    elif area == 'Leon_FL':
        x, y = base.get_xlim()[0]+75, base.get_ylim()[0]+75
        scale_text = '|     ~10km~     |'
    
    elif area == 'Test_Waverly_Leon_FL':
        x, y = base.get_xlim()[0]+75, base.get_ylim()[0]+75
        scale_text = '|    ~.5km~    |'
    
    elif area == 'Test_Grid_Leon_FL':
        x, y = base.get_xlim()[0]+15, base.get_ylim()[0]+15
        scale_text = '|    ~50m~    |'
    
    elif area == 'Test_Small_Leon_FL':
        x, y = base.get_xlim()[0]+15, base.get_ylim()[0]+15
        scale_text = '|    ~50m~    |'
    '''
    
    if for_multiplot:
        fontsize = 'small'
    else:
        fontsize = 'medium'
    scale_text = '|    ~%s%s~    |' % (distance, units)
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    base.text(x, y, scale_text, fontsize=fontsize,
              fontstyle='italic', bbox=bbox_props)


def add_legend(patches, for_multiplot=False):
    """Add a legend to a plot
    
    Parameters
    ----------
    patches : list
        legend handles matching plotted items
    for_multiplot : create_patches 
    """
    if for_multiplot:
        anchor = (1.1, 1.65)
    else:
        anchor = (1.005, 1.016)
    legend = plt.legend(handles=patches, loc='upper left',
                        fancybox=True, framealpha=.85,
                        bbox_to_anchor=anchor, fontsize='x-large')
    legend.get_frame().set_facecolor('white')


def dv_colorset(dvs):
    """decision variables color set
    
    Parameters
    ---------
    dvs : geopandas.GeoSeries
        facility decision variables
    
    Returns
    -------
    dv_colors : dict
        decision variable to set color lookup
    """
    dv_colors = ['fuchsia', 'mediumseagreen', 'blueviolet',
                 'darkslategray', 'lightskyblue', 'saddlebrown',
                 'cyan', 'darkgoldenrod', 'limegreen', 'peachpuff',
                 'coral', 'mediumvioletred', 'darkcyan',
                 'thistle', 'lavender', 'tomato']
    dv_colors = {desc_var:dv_colors[idx] for idx, desc_var\
                 in enumerate(dvs)}
    return dv_colors


def create_patches(model=None, pt1_size=None, pt2_size=None,
                   buffered=None, legend_aux=None, dvs_to_leg=None,
                   for_multiplot=False, style=None, census_geo=None):
    """create all patches to add to the legend.
    
    Parameters
    ----------
    for_multiplot : bool
        for a single plot (True), or multiplot (False).
        Default is False.
    model : see plotter()
    pt1_size : see plotter()
    pt2_size : see plotter()
    buffered : see plotter()
    census_geo : see ploter()
    legend_aux : see plotter()
    dvs_to_leg : see plotter()
    
    Returns
    -------
    patches : list
        legend handles matching plotted items
    """
    if pt1_size:
        ms1 = float(pt1_size)/6.
    if pt2_size:
        ms2 = float(pt2_size)/8.
    
    spacer = mpatches.Patch([], [], color='w', linewidth=0,
                            alpha=.0, label='')
    # all patches to add to legend
    patches = []
    # streets -- always plot
    street_width = 2
    strs = mlines.Line2D([], [], color='k', label='Streets',
                         alpha=1, linewidth=street_width)
    patches.extend([spacer, strs])
    # non-results data
    if legend_aux:
        if 'buffer' in legend_aux:
            label = 'Street buffer (%sm)' % buffered
            strbuff = mpatches.Patch([], [], color='y', linewidth=2,
                                     alpha=.5, label=label)
            patches.extend([spacer, strbuff])
        
        if 'census_polys' in legend_aux:
            label = 'Census %s ' % census_geo.capitalize()\
                    + '($n$=%s)' % model.n_cen
            cenpoly = mpatches.Patch([], [], color='b', alpha=.035,
                                     linestyle='--', linewidth=2,
                                     label=label)
            patches.extend([spacer, cenpoly])
        
        if 'cli_tru' in legend_aux:
            try:
                if dvs_to_leg:
                    pcolor = 'k'
                    msize = ms2/3.
                    plabel = 'Uncovered Households '\
                             + '($n$=%s)' % model.n_cli_uncov
                else:
                    pcolor = 'b'
                    msize = ms1
                    plabel = 'Households ($n$=%s)' % model.n_cli
                cli_tru = mlines.Line2D([], [], color=pcolor,
                                        marker='o', ms=msize,
                                        linewidth=0, alpha=1,
                                        markeredgecolor='k',
                                        label=plabel)
                patches.extend([spacer, cli_tru])
            except AttributeError:
                pass
        if 'fac_tru' in legend_aux:
            if dvs_to_leg:
                pcolor = 'k'
                msize = ms2
                pmarker = '*'
                no_fac = model.n_fac - len(list(model.fac2cli.keys()))
                plabel = 'Unselected Facilities ($n$=%s)' % no_fac
            else:
                pcolor = 'r'
                msize = ms1
                pmarker = 'o'
                plabel = 'Fire Stations'\
                         + '($n$=%s)' % model.n_fac
            fac_tru = mlines.Line2D([], [], color=pcolor,
                                    marker=pmarker, ms=msize,
                                    markeredgecolor='k',
                                    linewidth=0, alpha=1,
                                    label=plabel)
            patches.extend([spacer, fac_tru])
        if 'cli_snp' in legend_aux:
            label = 'Households snapped to network'
            cli_snp = mlines.Line2D([], [], color='b', marker='o',
                                    ms=ms2, linewidth=0, alpha=1,
                                    markeredgecolor='k', label=label)
            patches.extend([spacer, cli_snp])
        if 'fac_snp' in legend_aux:
            label = 'Fire Stations snapped to network'
            fac_snp = mlines.Line2D([], [], color='r', marker='o',
                                    ms=ms2, markeredgecolor='k',
                                    linewidth=0, alpha=1,
                                    label=label)   
            patches.extend([spacer, fac_snp])
    patches.extend([spacer])
    # results data for single plot
    if dvs_to_leg and not for_multiplot:
        # add facility, client, and service area patches to legend
        for k, v in dvs_to_leg.items():
            fdv_label = 'Fire Station %s' % k
            fdv = mlines.Line2D([], [], color=v['color'], marker='*',
                                ms=ms1/2., markeredgecolor='k',
                                linewidth=0, alpha=.8, label=fdv_label)
            cdv_label = 'Households served by %s ' % k \
                        + '($n$=%s)' % v['clients']
            cdv = mlines.Line2D([], [], color=v['color'], marker='o',
                                ms=ms1/6., markeredgecolor='k',
                                linewidth=0, alpha=.5, label=cdv_label)
            serv_label = '%s service area' % k
            serv = mpatches.Patch([], [], color=v['color'], linewidth=2,
                                  alpha=.25, label=serv_label)
            patches.extend([spacer, fdv, cdv, serv, spacer])
    # results data for multiplot
    if dvs_to_leg and for_multiplot:
        for idx, (k, v, n) in dvs_to_leg.items():
            fdv = mlines.Line2D([], [], color=v, marker='*', ms=ms1/2,
                                markeredgecolor='k', linewidth=0,
                                alpha=.8, label='%s ($n$=%s)' % (k,n))
            patches.extend([spacer, fdv, spacer])
    return patches


