"""hard coded plotting of dissertation figures
"""
import matplotlib
matplotlib.use('agg')


# standard library imports
import sys

# non-standard library imports
import geopandas as gpd
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib_scalebar.scalebar import ScaleBar

from shapely.geometry import Point, LineString, Polygon

# project library imports
from at1866_Master import spaghetti as spgh


def reject_outliers(data, m=3):
    """
    """
    bool_mask = abs(data - np.mean(data)) <= m * np.std(data)
    return data[bool_mask], bool_mask

spacer = mpatches.Patch([], [], color='w', linewidth=0,
                            alpha=.0, label='')


def ch2_lvd_pp2n(LVOR_PNT_RHO):
    """
    
    """
    
    # line Voronoi Diagram figure
    area_prefix = 'Test_Sine_'
    area_base = 'Leon_FL'
    area_suffix = '_2010.shp'
    area = area_prefix + area_base
    place_time_shp = area_base + area_suffix
    
    data = '%s/%s/intermediary/' % ('../data', area)
    net_dir = data
    
    data = '%s/%s/clean/' % ('../data', area)
    alc_dir = data + 'allocation_data/'
    
    segms = gpd.read_file('%sNetSegms_%s' % (net_dir, place_time_shp))
    nodes = gpd.read_file('%sNetNodes_%s' % (net_dir, place_time_shp))
    
    ss = segms[:11]
    sn = nodes[:8]
    
    lvd = gpd.read_file('%slvd_%s_Test_Sine_Leon_FL.shp' % (alc_dir, LVOR_PNT_RHO))
    
    lvd_handle = mpatches.Patch(fc='w', ec='r', alpha=.7, ls='-', lw=2.5, hatch='/',
                            label='Voronoi cells\n$\longrightarrow$ $R$=%s' % lvd.shape[0])
    street_handle = mlines.Line2D([], [], color='k', alpha=1, linewidth=3.5,
                                  label='Network segments\n$\longrightarrow$ $S$=%s' % ss.shape[0])
    vertex_handle = mlines.Line2D([], [], color='k', alpha=1, marker='o', ms=10,
                                  linewidth=0,
                                  label='Network vertices\n$\longrightarrow$ $V$=%s' % sn.shape[0])
    
    handles = [spacer, lvd_handle,
               spacer, street_handle,
               spacer, vertex_handle,
               spacer]
    
    fig, ax = plt.subplots(figsize=(9,9))
    
    lvd.plot(ax=ax, column='SegID', cmap='Blues', alpha=.7,
                    linestyle='-', edgecolor='r', linewidth=2.5)
    
    ss.plot(ax=ax, linewidth=3.5, color='k', alpha=1)
    
    sn.plot(ax=ax, markersize=100, color='k', alpha=1, zorder=2)
    
    plt.legend(handles=handles, bbox_to_anchor=(1., .875), fontsize=18);
    
    plt.axis('off')
    
    plt.savefig('../results/Test_Sine_Leon_FL/plots/line_voronoi.png',
                bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    ######################################################################
    
    # pp2n concept figure
    pop_block = [Polygon(((0,0), (4,0), (4,4), (0,4))),
                 Polygon(((4,0), (8,0), (8,4), (6,4), (6,1), (4,1), ))]
    pop_block_gdf = gpd.GeoDataFrame(geometry=pop_block)
    pop_block_handle = mpatches.Patch(fc='b', ec='k', alpha=.1,
                                      label='Populated geographies', ls='--', lw=2)
    
    nop_block = [Polygon(((4,1), (6,1), (6,4), (4,4)))]
    nop_block_gdf = gpd.GeoDataFrame(geometry=nop_block)
    nop_block_handle = mpatches.Patch(fc='w', ec='k', alpha=.1,
                                      label='Unpopulated geography', ls='--', lw=2)
    
    street = [LineString((Point(4,-1), Point(4,4.4)))]
    street_gdf = gpd.GeoDataFrame(geometry=street)
    street_handle = mlines.Line2D([], [], color='k', label='Network segment',
                                  alpha=1, linewidth=3)
    
    offsets = [LineString((Point(3.7,-1), Point(3.7,4.4))),
           LineString((Point(4.3,-1), Point(4.3,4.4)))]
    offsets_gdf = gpd.GeoDataFrame(geometry=offsets)
    offsets_handle = mlines.Line2D([], [], color='k',
                                   label='Step 2(c) parallel offset',
                                   alpha=.5, ls=':', linewidth=2)
    
    points_good = [Point(3.7,2), Point(4.3,0.8)]
    points_good_gdf = gpd.GeoDataFrame(geometry=points_good)
    points_good_handle = mlines.Line2D([], [], color='g', label='Step 2(d) successes',
                                        alpha=1, marker='o', ms=10,
                                        linewidth=0)
    
    points_bad = [Point(4.3,2), Point(4.3,1.4),
                  Point(4.3,2.6), Point(4.3,3.2),]
    points_bad_gdf = gpd.GeoDataFrame(geometry=points_bad)
    points_bad_handle = mlines.Line2D([], [], color='r', label='Step 2(d) failures',
                                        alpha=1, marker='x', ms=12.5,
                                        linewidth=0)
    
    handles = [pop_block_handle, spacer,
           nop_block_handle, spacer,
           street_handle, spacer,
           offsets_handle, spacer,
           points_good_handle, spacer,
           points_bad_handle]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    pop_block_gdf.plot(ax=ax, color='b', alpha=.1, edgecolor='k',
                       linestyle='--', lw=2)
    
    nop_block_gdf.plot(ax=ax, color='w', alpha=.1, edgecolor='k',
                       linestyle='--', lw=2)
    
    street_gdf.plot(ax=ax, color='k', alpha=1, edgecolor='k', lw=3)
    
    offsets_gdf.plot(ax=ax, color='k', alpha=.5, linestyle=':', edgecolor='k', lw=3)
    
    points_good_gdf.plot(ax=ax, color='g', alpha=1, markersize=100, zorder=2)
    
    points_bad_gdf.plot(ax=ax, color='r', alpha=1, markersize=200, zorder=2, marker='x')
    
    ax.annotate('1', xy=(0,0), xytext=(3.3, 1.9), size=15)
    ax.annotate('1', xy=(0,0), xytext=(4.6, 1.9), size=15)
    ax.annotate('2', xy=(0,0), xytext=(4.6, 2.5), size=15)
    ax.annotate('3', xy=(0,0), xytext=(4.6, 1.3), size=15)
    ax.annotate('4', xy=(0,0), xytext=(4.6, 3.1), size=15)
    ax.annotate('5', xy=(0,0), xytext=(4.6, 0.7), size=15)
    
    ax.annotate('$s_{i}$', xy=(0,0), xytext=(3.9, 4.45), size=17)
    ax.annotate('$s_{iL}$', xy=(0,0), xytext=(3.3, 4.2), size=17)
    ax.annotate('$s_{iR}$', xy=(0,0), xytext=(4.4, 4.2), size=17)
    
    ax.annotate('$pnt_{i}$', xy=(0,0), xytext=(2.8, 1.6), size=15)
    ax.annotate('$pnt_{i+1}$', xy=(0,0), xytext=(4.6, .4), size=15)
    
    plt.legend(handles=handles, bbox_to_anchor=(1.65, .97), fontsize=18)
    
    plt.axis('off')
    
    plt.savefig('../results/Test_Sine_Leon_FL/plots/pp2n_concept.png',
                bbox_inches='tight', edgecolor='k',
                format='png', dpi=400 ,quality=95)


def ch2_grid_plots(adjusted=False):
    """
    Chapter 2 grid 1 - 
    Simulated data on a Cartesian plane. Household information includes 
    true location, network proximity, and total residents per 
    household. Segments are labeled.
    
    Chapter 2 grid 2 - 
    Derived population representations. The \textit{pc2n} and 
    \textit{pp2n} representations, as well as the true locations must be 
    snapped to network segments, whereas the \textit{va2n} representations are
     generated on network segments. All methods for representing population 
     are valid in that they sum to the total true population of 48.
    
    """
    if adjusted:
        tag = '_ADJUSTED'
    else:
        tag = ''
    
    ################################## grid 1
    area_prefix = 'Test_Grid_'
    area_base = 'Leon_FL'
    area_suffix = '_2010' + tag + '.shp'
    area = area_prefix + area_base
    place_time_shp = area_base + area_suffix
    
    
    data = '%s/%s/clean/' % ('../data', area)
    
    obs_dir = data + 'observation_data/'
    cen_dir = data + 'census_data/'
    net_dir = data + 'network_data/'
    mtx_dir = data + 'cost_matrices/'
    alc_dir = data + 'allocation_data/'
    
    blks_file = '%sCensusBlocks_%s' % (cen_dir, place_time_shp)
    blks = gpd.read_file(blks_file)
    blocks_pop = blks[blks.synth_pop > 0]
    blocks_nop = blks[blks.synth_pop == 0]

    pop_block_handle = mpatches.Patch(fc='b', ec='k', alpha=.1,
                                      label='Populated geographies', ls='--', lw=2)
    nop_block_handle = mpatches.Patch(fc='w', ec='k', alpha=.1,
                                      label='Unpopulated geographies', ls='--', lw=2)

    # households
    hh_file = '%sHouseholds_Synthetic_%s' % (cen_dir, place_time_shp)
    hh_gdf = gpd.read_file(hh_file)
    hh_handle = mlines.Line2D([], [], color='w', markeredgecolor='k',
                              label='Weight-labeled\nsynthetic households',
                                            alpha=1, marker='s', ms=12,
                                            linewidth=0)

    streets_file = '%sSimplifiedSegms_%s' % (net_dir, place_time_shp)
    streets = gpd.read_file(streets_file)
    street_handle = mlines.Line2D([], [], color='k',
                                  marker='o', markersize=10, 
                                  markeredgecolor='k',
                                  markerfacecolor='PALEGOLDENROD',
                                  label='ID-labeled\nNetwork segments',
                                     alpha=1., linewidth=2)

    vertexx_file = '%sSimplifiedNodes_%s' % (net_dir, place_time_shp)
    vertexx = gpd.read_file(vertexx_file)
    vertex_handle = mlines.Line2D([], [], color='k', label='Network vertices',
                                            alpha=1, marker='o', ms=5,
                                            linewidth=0)

    handles = [spacer,pop_block_handle,
               spacer,nop_block_handle,
               spacer,hh_handle,
               spacer,street_handle,
               spacer,vertex_handle,
               spacer]

    fig, ax = plt.subplots(figsize=(9, 9))
    blocks_pop.plot(ax=ax, color='b', alpha=.1, edgecolor='k',
                       linestyle='--', lw=2)
    blocks_nop.plot(ax=ax, color='w', alpha=.1, edgecolor='k',
                       linestyle='--', lw=2)

    hh_gdf.plot(ax=ax, color='w', edgecolor='k', alpha=1, marker='s', markersize=200, zorder=2)

    streets.plot(ax=ax, color='k', alpha=1., lw=2)

    vertexx.plot(ax=ax, color='k', alpha=1, markersize=20, zorder=2)

    # weights
    ax.annotate('12', xy=(0,0), xytext=(6.03, 8.7555), size=11)
    ax.annotate('9', xy=(0,0), xytext=(7.465, 4.485), size=11)
    ax.annotate('9', xy=(0,0), xytext=(8.72, 1.85), size=11)
    ax.annotate('6', xy=(0,0), xytext=(5.12, 7.55), size=11)
    ax.annotate('6', xy=(0,0), xytext=(3.75, 5.17), size=11)
    ax.annotate('3', xy=(0,0), xytext=(4.41, 2.625), size=11)
    ax.annotate('3', xy=(0,0), xytext=(.175, 8.68), size=11)

    # segments
    seg_font_size = 11
    ax.annotate('0', xy=(0,0), xytext=(2.92, 1.425), size=seg_font_size)
    ax.annotate('1', xy=(0,0), xytext=(2.92, 4.425), size=seg_font_size)
    ax.annotate('2', xy=(0,0), xytext=(2.92, 7.425), size=seg_font_size)
    ax.annotate('3', xy=(0,0), xytext=(5.92, 1.425), size=seg_font_size)
    ax.annotate('4', xy=(0,0), xytext=(5.92, 4.425), size=seg_font_size)
    ax.annotate('5', xy=(0,0), xytext=(5.92, 7.425), size=seg_font_size)

    ax.annotate('6', xy=(0,0), xytext=(1.425, 2.925), size=seg_font_size)
    ax.annotate('7', xy=(0,0), xytext=(1.425, 5.925), size=seg_font_size)

    ax.annotate('8', xy=(0,0), xytext=(4.425, 2.925), size=seg_font_size)
    ax.annotate('9', xy=(0,0), xytext=(4.425, 5.925), size=seg_font_size)

    ax.annotate('10', xy=(0,0), xytext=(7.35, 2.9), size=seg_font_size)
    ax.annotate('11', xy=(0,0), xytext=(7.35, 5.9), size=seg_font_size)

    to_buffer = .15
    circles = [Point(3, 1.5).buffer(to_buffer),#0
               Point(3, 4.5).buffer(to_buffer),#1
               Point(3, 7.5).buffer(to_buffer),#2
               Point(6, 1.5).buffer(to_buffer),#3
               Point(6, 4.5).buffer(to_buffer),#4
               Point(6, 7.5).buffer(to_buffer),#5
               Point(1.5, 3).buffer(to_buffer),#6
               Point(1.5, 6).buffer(to_buffer),#7
               Point(4.5, 3).buffer(to_buffer),#8
               Point(4.5, 6).buffer(to_buffer),#9
               Point(7.5, 3).buffer(to_buffer),#10
               Point(7.5, 6).buffer(to_buffer) #11
               ]

    circles_gdf = gpd.GeoDataFrame(geometry=circles)
    circles_gdf.plot(ax=ax, edgecolor='k', facecolor='PALEGOLDENROD', zorder=3)

    y_factor = .9
    plt.legend(handles=handles, bbox_to_anchor=(1.65, y_factor), fontsize=18)

    plt.savefig('../results/Test_Grid_Leon_FL/plots/grid_1'+tag+'.png',
                bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    plt.close()


    ###################################################### grid 2
    pp = 'PointPattern'
    census_bases = ['pc2nPopulatedBlocks'+tag,
                    'HouseholdsSynthetic']
    allocation_bases = ['va2nPopulatedBlocks'+tag,
                        'pp2n0.25PopulatedBlocks'+tag]

    patterns = {name: spgh.load_pickled(cen_dir,
                                        pickle_name='%s_%s' % (pp, name))\
                                        for name in census_bases}

    patterns.update({name: spgh.load_pickled(alc_dir,
                                             pickle_name='%s_%s' % (pp, name))\
                                             for name in allocation_bases})
                                         
    pc2n = patterns['pc2nPopulatedBlocks'+tag].df
    pp2n = patterns['pp2n0.25PopulatedBlocks'+tag].df
    va2n = patterns['va2nPopulatedBlocks'+tag].df

    pc2n_handle = mlines.Line2D([], [], color='g', marker='o',
                                label='$pc2n$\n$\longrightarrow$ $RP$=%s, $\sum pop$=%s' % (pc2n.shape[0], pc2n.synth_pop.sum()),
                                alpha=1, ms=7, linewidth=0)

    va2n_handle = mlines.Line2D([], [], color='r', marker='o',
                                label='$va2n$\n$\longrightarrow$ $RP$=%s, $\sum pop$=%s' % (va2n.shape[0], round(va2n.pop_va2n.sum(), 3)),
                                alpha=1, ms=7, linewidth=0)
    
    pp2n_handle = mlines.Line2D([], [], color='b', marker='o',
                                label='$pp2n$ (0.25)\n$\longrightarrow$ $RP$=%s, $\sum pop$=%s' % (pp2n.shape[0], pp2n.pop_pp2n.sum()),
                                alpha=1, ms=7, linewidth=0)

    hh_handle = mlines.Line2D([], [], color='w', markeredgecolor='k',
                              label='Synthetic households\n$\longrightarrow$ $RP$=%s, $\sum pop$=%s' % (hh_gdf.shape[0], hh_gdf.synth_pop.sum()),
                                            alpha=1, marker='s', ms=12,
                                            linewidth=0)

    handles = [pop_block_handle,
               spacer,nop_block_handle,
               spacer,hh_handle,
               spacer,pc2n_handle,
               spacer,va2n_handle,
               spacer,pp2n_handle,
               spacer,street_handle]
    
    fig, ax = plt.subplots(figsize=(9, 9))

    blocks_pop.plot(ax=ax, color='b', alpha=.1, edgecolor='k',
                       linestyle='--', lw=2)
    blocks_nop.plot(ax=ax, color='w', alpha=.1, edgecolor='k',
                       linestyle='--', lw=2)

    streets.plot(ax=ax, color='k', alpha=.5, lw=2)

    hh_gdf.plot(ax=ax, color='w', edgecolor='k', alpha=1, marker='s', markersize=200, zorder=2)

    pc2n.plot(ax=ax, color='g', alpha=1, markersize=50, zorder=2)
    pp2n.plot(ax=ax, color='b', alpha=1, markersize=50, zorder=2)
    va2n.plot(ax=ax, color='r', alpha=1, markersize=50, zorder=2)
    
    y_factor = .97
    plt.legend(handles=handles, bbox_to_anchor=(1.65, y_factor), fontsize=18)

    plt.savefig('../results/Test_Grid_Leon_FL/plots/grid_2'+tag+'.png',
                bbox_inches='tight',
                format='png', dpi=400 ,quality=95)

    plt.close()



def ch2_tract_plots(adjusted=False, to_fs=False):
    """
    Chapter 2 tract 1 - 
    2010 Decennial census geographies in Tract 12073001700 of Leon County, FL.
    
    Chapter 2 tract 2 - 
    2011(2010-certified) population-weighted property parcel 
        centroids within Tract 12073001700 in Leon County, FL and the 
        intersecting network structure.
    
    Chapter 2 tract zoom -
    The flexibility and robustness of the \textit{pp2n} method is 
         demonstrated by restricting the generation of points to only 
         segments which are associated with populated census geographies 
         and within the defined offset parameter.
    
    Chapter 2 euclidean -- network -
    Network distance and the analogous ratio to Euclidean distance 
         for unique, randomly sampled origin-destination (OD) pairs 
         between weighted property parcels and abstract population 
         representations generated from populated census blocks. 
         Each plot axis is annotated with $N$ total observation 
         and $n$ sample observations.
    """
    
    if adjusted:
        tag = '_ADJUSTED'
    else:
        tag = ''
    
    area_prefix = 'Test_Tract_'
    area_base = 'Leon_FL'
    area_suffix = '_2010'+tag+'.shp'
    area = area_prefix + area_base
    place_time_shp = area_base + area_suffix
    
    data = '%s/%s/clean/' % ('../data', area)
    
    obs_dir = data + 'observation_data/'
    cen_dir = data + 'census_data/'
    net_dir = data + 'network_data/'
    mtx_dir = data + 'cost_matrices/'
    alc_dir = data + 'allocation_data/'
    
    
    ############################################### county context
    county_file = '../%s/%s/%s/%s/Counties_%s' % ('data', area_base, 'clean', 'census_data',place_time_shp)
    county = gpd.read_file(county_file)
    
    city_file = '%sTallahassee_%s' % (cen_dir, place_time_shp)
    city = gpd.read_file(city_file)
    
    all_tracts_file = '../%s/%s/%s/%s/CensusTracts_%s' % ('data', area_base, 'clean', 'census_data',place_time_shp)
    all_tracts = gpd.read_file(all_tracts_file)
    
    tracts_file = '%sCensusTracts_%s' % (cen_dir, place_time_shp)
    tracts = gpd.read_file(tracts_file)
    
    one_tract_lbl = mpatches.Patch(fc='gray', ec='k', alpha=1., lw=3,
                                    label='Tract 12073001700')
    trtlbl = mpatches.Patch(fc='w', ec='k', alpha=.2, lw=1,
                                label='Tracts\n$\longrightarrow$ $POLY$=%s' % all_tracts.shape[0])
    city_handle = mpatches.Patch(fc='gold', ec='gold',
                                  label='Tallahassee city limits',
                                     alpha=.5, linewidth=1)
    county_handle = mlines.Line2D([], [], color='k',
                                      label='Leon County boundary',
                                         alpha=.75, linewidth=1)
    
    handles = [one_tract_lbl, trtlbl, county_handle, city_handle]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    county.plot(ax=ax, color='w', alpha=1, lw=1, edgecolor='k', zorder=0)
    city.plot(ax=ax, color='gold', alpha=.5, lw=1, edgecolor='gold', zorder=0)
    all_tracts.plot(ax=ax, color='w', alpha=.2, lw=1, edgecolor='k', zorder=0)
    tracts.plot(ax=ax, color='gray',alpha=1., edgecolor='k', lw=3, zorder=1)
    
    plt.legend(handles=handles,
                   bbox_to_anchor=(0., -0., 1., -0.),
                   fontsize=15,
                   ncol=2, mode="expand",
                   borderpad=1.,
                   columnspacing=1.,
                   labelspacing=1.)
    
    x, y = 651000, 144000
    arw = 'rarrow, pad=0.15'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.15, location='lower left')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 578000, 186000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Test_Tract_Leon_FL/plots/tract_leon_full_ch2.png', bbox_inches='tight',
                    format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    ################################################## tract 1
    area_prefix = 'Test_Tract_'
    area_base = 'Leon_FL'
    area_suffix = '_2010'+tag+'.shp'
    area = area_prefix + area_base
    place_time_shp = area_base + area_suffix
    
    data = '%s/%s/clean/' % ('../data', area)
    
    obs_dir = data + 'observation_data/'
    cen_dir = data + 'census_data/'
    net_dir = data + 'network_data/'
    mtx_dir = data + 'cost_matrices/'
    alc_dir = data + 'allocation_data/'
    
    blks_file = '%sCensusBlocks_%s' % (cen_dir, place_time_shp)
    blks = gpd.read_file(blks_file)
    blocks_pop = blks[blks.POP100 > 0]
    blocks_nop = blks[blks.POP100 == 0]
    
    tracts_file = '%sCensusTracts_%s' % (cen_dir, place_time_shp)
    tracts = gpd.read_file(tracts_file)
    bkgs_file = '%sCensusBlockGroups_%s' % (cen_dir, place_time_shp)
    bkgs = gpd.read_file(bkgs_file)
    
    trtlbl = mpatches.Patch(fc='w', ec='k', alpha=.2, lw=5,
                            label='Tracts\n$\longrightarrow$ $POLY$=%s' % tracts.shape[0])
    bgrlbl = mpatches.Patch(fc='w', ec='k', alpha=.3, lw=2,
                            label='Block Groups\n$\longrightarrow$ $POLY$=%s' % bkgs.shape[0])
    pbklbl = mpatches.Patch(fc='w', ec='k', alpha=.2, ls='-.',
                            label='Populated Blocks\n$\longrightarrow$ $POLY$=%s' % blocks_pop.shape[0])
    nbklbl = mpatches.Patch(fc='r', ec='k', alpha=.75, ls='-.',
                            label='Unpopulated Blocks\n$\longrightarrow$ $POLY$=%s' % blocks_nop.shape[0])
    
    handles = [spacer,trtlbl,
               spacer,bgrlbl,
               spacer,pbklbl,
               spacer,nbklbl,spacer]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    tracts.plot(ax=ax, color='w',alpha=.2, edgecolor='k', lw=10)
    bkgs.plot(ax=ax, color='w', alpha=.1, edgecolor='k', lw=4)
    blocks_pop.plot(ax=ax, color='w', alpha=.2, edgecolor='k',  linestyle='-.')
    blocks_nop.plot(ax=ax, color='r', alpha=.75, edgecolor='k',  linestyle='-.')
    
    plt.legend(handles=handles, bbox_to_anchor=(1.6, .85), fontsize=18);
    
    ax.set_xlim(620800, 625600)
    ax.set_ylim(161900, 166800)
    
    x, y = 625250, 166200
    arw = 'rarrow, pad=0.25'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.3, location='lower right')
    ax.add_artist(scalebar)

    bbox_props = dict(boxstyle='round, pad=.5',
                  fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 620900, 162000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Test_Tract_Leon_FL/plots/tract_geographies_ch2'+tag+'.png', bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    plt.close()
    
    
    ##########################################################################
    ############### tract 2
    # parcels centroids
    parcel_centroids_file = '%sWeightedParcels_%s' % (cen_dir,place_time_shp)
    parcel_centroids = gpd.read_file(parcel_centroids_file)
    
    pop_parcel_centroids = parcel_centroids[parcel_centroids['SUM_EST_PO'] > 0.0]
    no_pop_parcel_centroids = parcel_centroids[parcel_centroids['SUM_EST_PO'] <= 0.]
    
    pop_hdl = mlines.Line2D([], [], color='m', markeredgecolor='k',
                        alpha=.5, marker='s', ms=5, linewidth=0,
                            label='Weighted Parcels\n$\\rightarrow$ $RP$=%s' % pop_parcel_centroids.shape[0]\
                                + '\n$\\rightarrow$ $\sum pop$=%s' % round(pop_parcel_centroids['SUM_EST_PO'].sum(), 8))
    
    no_pop_hdl = mlines.Line2D([], [], color='r', markeredgecolor='k',
                            alpha=1., marker='s', ms=6.5, linewidth=0,
                                label='Unweighted Parcels\n$\\rightarrow$ $RP$=%s' % no_pop_parcel_centroids.shape[0]\
                                       + '\n$\\rightarrow$ $\sum pop$=%s' % no_pop_parcel_centroids['SUM_EST_PO'].sum())
    
    streets_file = '%sSimplifiedSegms_%s' % (net_dir, place_time_shp)
    streets = gpd.read_file(streets_file)
    street_handle = mlines.Line2D([], [], color='k',
                                  label='Network segments\n$\longrightarrow$ $S$=%s' % streets.shape[0],
                         alpha=.7, linewidth=1)
    
    vertexx_file = '%sSimplifiedNodes_%s' % (net_dir, place_time_shp)
    vertexx = gpd.read_file(vertexx_file)
    vertex_handle = mlines.Line2D([], [], color='k',
                                  label='Network vertices\n$\longrightarrow$ $V$=%s' % vertexx.shape[0],
                                            alpha=1, marker='o', ms=2,
                                            linewidth=0)
    
    handles = [
           spacer,pop_hdl,
           spacer,no_pop_hdl,
           spacer,street_handle,
           spacer, vertex_handle, spacer]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    pop_parcel_centroids.plot(ax=ax, color='m', edgecolor='k',
                       alpha=.5, marker='s', markersize=6.5, zorder=2)
    no_pop_parcel_centroids.plot(ax=ax, color='r', edgecolor='k', lw=.25,
                       alpha=1, marker='s', markersize=10.5, zorder=2)
    streets.plot(ax=ax, color='k', alpha=.5, lw=1)
    vertexx.plot(ax=ax, color='k', alpha=1, markersize=2.5, zorder=1)
    
    plt.legend(handles=handles, bbox_to_anchor=(1.05, .9), fontsize=18)
    
    ax.set_xlim(620800, 625600)
    ax.set_ylim(161900, 166800)
    
    x, y = 625250, 166200
    arw = 'rarrow, pad=0.25'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    
    '''
    distance, units = 1000, 'm'
    scale_text = '|     ~  %s %s  ~    |' % (distance, units)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    
    ax.text(623975, 162100, scale_text, fontsize=10,
              fontstyle='italic', bbox=bbox_props)
    '''
    
    #######################################################################################################
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.3, location='lower right')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                  fc='w', ec='0.5', alpha=0.7)
    #############################################################################################################
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 620900, 162000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Test_Tract_Leon_FL/plots/tract_parcel_centroids_ch2'+tag+'.png', bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    plt.close()
    
    
    ##########################################################################
    ############### tract 2
    # parcels 
    
    parcel_file = '%sParcelPolygons_%s' % (cen_dir,place_time_shp)
    tract_parcels = gpd.read_file(parcel_file)
    
    no_pop_parcels = tract_parcels[tract_parcels['SUM_EST_PO'] <= 0.]
    pop_parcels = tract_parcels[tract_parcels['SUM_EST_PO'] > 0.0]
    
    # ambiguous population
    ambig_pop_parcels = no_pop_parcels[no_pop_parcels.PARCEL_ID\
                               .isin(no_pop_parcel_centroids.PARCEL_ID)]
    
    # absolutely no populated
    absol_no_pop_parcels = no_pop_parcels[~no_pop_parcels.PARCEL_ID\
                                      .isin(no_pop_parcel_centroids.PARCEL_ID)]
    
    
    pop_hdl = mpatches.Patch(fc='w', ec='k', alpha=1., lw=.1,
                            label='Populated Parcels\n$\\rightarrow$ $POLY$=%s' % pop_parcels.shape[0]\
                                + '\n$\\rightarrow$ $\sum pop$=%s' % round(pop_parcels['SUM_EST_PO'].sum(), 8))
    
    ambig_hdl = mpatches.Patch(fc='w', ec='g', alpha=.5, lw=.1, hatch='//////',
                            label='Ambiguous Parcels\n$\\rightarrow$ $POLY$=%s' % ambig_pop_parcels.shape[0]\
                                + '\n$\\rightarrow$ $\sum pop$=%s' % ambig_pop_parcels['SUM_EST_PO'].sum())
    
    no_pop_hdl = mpatches.Patch(fc='r', ec='k', alpha=.75, lw=.1,
                            label='Unpopulated Parcels\n$\\rightarrow$ $POLY$=%s' % absol_no_pop_parcels.shape[0]\
                                   + '\n$\\rightarrow$ $\sum pop$=%s' % round(absol_no_pop_parcels['SUM_EST_PO'].sum(),20))
    
    handles=[spacer,
         pop_hdl, spacer,
         ambig_hdl, spacer,
         no_pop_hdl, spacer]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    pop_parcels.plot(ax=ax, color='w', alpha=1, lw=.25, edgecolor='k', zorder=0)
    ambig_pop_parcels.plot(ax=ax, color='w', alpha=.5, lw=.1, hatch='//////',
                        edgecolor='g', zorder=1)
    absol_no_pop_parcels.plot(ax=ax, color='r', alpha=.75, lw=.25,
                        edgecolor='k', zorder=1)
    
    plt.legend(handles=handles, bbox_to_anchor=(1.6, .835), fontsize=18)
    
    ax.set_xlim(620800, 625600)
    ax.set_ylim(161900, 166800)
    
    x, y = 625250, 166200
    arw = 'rarrow, pad=0.25'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.3, location='lower right')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                  fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 620900, 162000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Test_Tract_Leon_FL/plots/tract_parcels_ch2'+tag+'.png', bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    
    plt.close()



    #----------------------------------------
    # zoom
    pp = 'PointPattern'
    census_bases = ['pc2nPopulatedBlocks'+tag,
                    'WeightedParcels']
    allocation_bases = ['va2nPopulatedBlocks'+tag,
                        'pp2n5.0PopulatedBlocks'+tag,
                        'pp2n36.0PopulatedBlocks'+tag]
    
    patterns = {name: spgh.load_pickled(cen_dir,
                                        pickle_name='%s_%s' % (pp, name))\
                                        for name in census_bases}
    
    patterns.update({name: spgh.load_pickled(alc_dir,
                                             pickle_name='%s_%s' % (pp, name))\
                                             for name in allocation_bases})
    
    pc2n = patterns['pc2nPopulatedBlocks'+tag].df
    pp2n5 = patterns['pp2n5.0PopulatedBlocks'+tag].df
    pp2n36 = patterns['pp2n36.0PopulatedBlocks'+tag].df
    va2n = patterns['va2nPopulatedBlocks'+tag].df
    
    pc2n_handle = mlines.Line2D([], [], color='g', marker='o',
                                label='$pc2n$',
                                alpha=1, ms=7, linewidth=0)
    
    pp2n5_handle = mlines.Line2D([], [], color='b', marker='o',
                                label='$pp2n$ (5m)',
                                alpha=1, ms=7, linewidth=0)
    
    pp2n36_handle = mlines.Line2D([], [], color='b', marker='^',
                                label='$pp2n$ (36m)',
                                alpha=1, ms=7, linewidth=0)
    
    va2n_handle = mlines.Line2D([], [], color='r', marker='o',
                                label='$va2n$',
                                alpha=1, ms=7, linewidth=0)
    
    streets_file = '%sSimplifiedSegms_%s' % (net_dir, place_time_shp)
    streets = gpd.read_file(streets_file)
    street_handle = mlines.Line2D([], [], color='k',
                                  label='Network segments',
                         alpha=.7, linewidth=1)
    
    vertexx_file = '%sSimplifiedNodes_%s' % (net_dir, place_time_shp)
    vertexx = gpd.read_file(vertexx_file)
    vertex_handle = mlines.Line2D([], [], color='k',
                                  label='Network vertices',
                                            alpha=1, marker='o', ms=2,
                                            linewidth=0)
    
    parcel_handle = mlines.Line2D([], [], color='m', markeredgecolor='k',
                              label='$wp2n$',
                                alpha=.5, marker='s', ms=5,
                                            linewidth=0)
    
    pbklbl = mpatches.Patch(fc='k', ec='k', alpha=.05, ls='-.',
                            label='Populated Blocks')
    
    handles = [pbklbl, spacer,
               pc2n_handle, spacer,
               va2n_handle, spacer,
               pp2n5_handle, spacer,
               pp2n36_handle, spacer,
               parcel_handle, spacer,
               street_handle, spacer,
               vertex_handle]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    
    blocks_pop.plot(ax=ax, color='k', alpha=.05, edgecolor='k',  linestyle='-.')
    parcel_centroids.plot(ax=ax, color='m', edgecolor='k',
                       alpha=.5, marker='s', markersize=12, zorder=2)
    pc2n.plot(ax=ax, color='g', alpha=1, markersize=12, zorder=2)
    pp2n5.plot(ax=ax, color='b', alpha=1, markersize=12, zorder=2)
    pp2n36.plot(ax=ax, color='b', marker='^', alpha=1, markersize=12, zorder=2)
    va2n.plot(ax=ax, color='r', alpha=1, markersize=12, zorder=2)
    
    streets.plot(ax=ax, color='k', alpha=.5, lw=1, zorder=0)
    vertexx.plot(ax=ax, color='k', alpha=1, markersize=2.5, zorder=1)
    
    plt.legend(handles=handles, bbox_to_anchor=(1.045, -0.025), fontsize=15)
    
    ax.set_xlim(623200, 623800)
    ax.set_ylim(166000, 166400)
    
    x, y = 623250, 166325
    arw = 'rarrow, pad=0.25'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.2, location='lower right')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                  fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 623210, 166384
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Test_Tract_Leon_FL/plots/tract_zoom_ch2'+tag+'.png', bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    #----------------------------------------
    # orphan
    
    orphan_df = pp2n36[pp2n36['GEOID'] == '120730017002023']
    orphan_handle = mlines.Line2D([], [], color='k', marker='^',
                            label='$pp2n$ (36m) orphan',
                            alpha=1, ms=15, linewidth=0)
    
    handles = [pbklbl, spacer,
               pp2n5_handle, spacer,
               pp2n36_handle, spacer,
               orphan_handle, spacer,
               parcel_handle, spacer,
               street_handle, spacer,
               vertex_handle]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    
    blocks_pop.plot(ax=ax, color='k', alpha=.05, edgecolor='k',  linestyle='-.')
    parcel_centroids.plot(ax=ax, color='m', edgecolor='k',
                           alpha=.5, marker='s', markersize=12, zorder=2)
    pp2n5.plot(ax=ax, color='b', alpha=1, markersize=12, zorder=2)
    pp2n36.plot(ax=ax, color='b', marker='^', alpha=1, markersize=12, zorder=2)
    orphan_df.plot(ax=ax, color='k', marker='^', alpha=1, markersize=50, zorder=3)
    streets.plot(ax=ax, color='k', alpha=.5, lw=1, zorder=0)
    vertexx.plot(ax=ax, color='k', alpha=1, markersize=2.5, zorder=1)
    
    plt.legend(handles=handles, bbox_to_anchor=(1.425, -0.035), fontsize=15)
    
    plt.locator_params(axis='x', nbins=6)
    ax.set_xlim(623550, 623775)
    ax.set_ylim(164730, 164860)
    
    x, y = 623560, 164750
    arw = 'rarrow, pad=0.25'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.3, location='lower right')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                  fc='w', ec='0.5', alpha=0.7)
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 623712, 164853
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Test_Tract_Leon_FL/plots/tract_orphan_ch2'+tag+'.png', bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    #----------------------------------------
    # euclidean - network
    method1 = 'euclidean'
    method2 = 'network'
    unit1 = 'pc2nPopulatedBlocks'
    unit2 = 'va2nPopulatedBlocks'
    unit3 = 'WeightedParcels'
    unit4 = 'pp2n5.0PopulatedBlocks'
    unit5 = 'pp2n36.0PopulatedBlocks'
    to = '_to_'
    cmtx = '_DataFrame.csv'
    
    units = [unit1, unit2, unit3, unit4, unit5]
    unit_dict = {}
    for unit in units:
        
        if to_fs:
            unit_ = 'FireStationsSynthetic'
        else:
            unit_ = unit
        
        # file names
        euc_file = method1 + '_' + unit_  + to + unit + cmtx
        net_file = method2 + '_' + unit_  + to + unit + cmtx
        
        # read files
        net_matrix = pd.read_csv(mtx_dir+net_file, index_col=0)
        net_matrix = net_matrix.values
        
        euc_matrix = pd.read_csv(mtx_dir+euc_file, index_col=0)
        euc_matrix = euc_matrix.values
        
        # calc ratio
        ratio_matrix = net_matrix / euc_matrix
        # create vectors
        tri_nidx = net_matrix.shape[0]
        net_vector = net_matrix[np.triu_indices(tri_nidx, k=1)]
        rat_vector = ratio_matrix[np.triu_indices(tri_nidx, k=1)]
        
        total_size = rat_vector.shape[0]
        
        # remove outliers more than 1 std deviation away from mean
        rat_vector, mask = reject_outliers(rat_vector)
        net_vector = net_vector[mask]
        
        std3 = net_vector.shape[0]
        
        # add to dict
        unit_dict[unit] = {'net': net_vector,
                            'rat': rat_vector,
                            'OD': total_size,
                            'std3': std3}
                       
    
    gridsize = 25
    
    fig, arr = plt.subplots(3,2,figsize=(15,15),
                            sharex=True, sharey=True)
    
    # pc2n
    ax = arr[0,0]
    hb = ax.hexbin(unit_dict[unit1]['net'],
                   unit_dict[unit1]['rat'],
                    gridsize=gridsize, cmap='Greens', mincnt=1)
    ax.set_title('$pc2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # va2n
    ax = arr[0,1]
    hb = ax.hexbin(unit_dict[unit2]['net'],
                  unit_dict[unit2]['rat'], 
                    gridsize=gridsize, cmap='Reds', mincnt=1)
    ax.set_title('$va2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # parcels
    ax = arr[1,0]
    hb = ax.hexbin(unit_dict[unit3]['net'],
                  unit_dict[unit3]['rat'], 
                    gridsize=gridsize, cmap='RdPu', mincnt=1)
    ax.set_title('$wp2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    
    if to_fs:
        arr[1,1].text(500, 4.35, 
                      '$fs$ $\\rightarrow$ $pc2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit1]['std3'], unit_dict[unit1]['OD']), fontsize=18)
        arr[1,1].text(500, 3.5, 
                      '$fs$ $\\rightarrow$ $va2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit2]['std3'], unit_dict[unit2]['OD']), fontsize=18)
        arr[1,1].text(500, 2.65, 
                      '$fs$ $\\rightarrow$ $wp2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit3]['std3'], unit_dict[unit3]['OD']), fontsize=18)
        arr[1,1].text(500, 1.8, 
                      '$fs$ $\\rightarrow$ $pp2n$ (5m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit4]['std3'], unit_dict[unit4]['OD']), fontsize=18)
        arr[1,1].text(500, 0.95, 
                      '$fs$ $\\rightarrow$ $pp2n$ (36m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit5]['std3'], unit_dict[unit5]['OD']), fontsize=18)
    else:
        arr[1,1].text(500, 4.35, 
                      '$pc2n$ $\\rightarrow$ $pc2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit1]['std3'], unit_dict[unit1]['OD']), fontsize=18)
        arr[1,1].text(500, 3.5, 
                      '$va2n$ $\\rightarrow$ $va2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit2]['std3'], unit_dict[unit2]['OD']), fontsize=18)
        arr[1,1].text(500, 2.65, 
                      '$wp2n$ $\\rightarrow$ $wp2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit3]['std3'], unit_dict[unit3]['OD']), fontsize=18)
        arr[1,1].text(500, 1.8, 
                      '$pp2n$ $\\rightarrow$ $pp2n$ (5m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit4]['std3'], unit_dict[unit4]['OD']), fontsize=18)
        arr[1,1].text(500, 0.95, 
                      '$pp2n$ $\\rightarrow$ $pp2n$ (36m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit5]['std3'], unit_dict[unit5]['OD']), fontsize=18)
    
    arr[1,1].axis('off')
    
    # pp2n 5
    ax = arr[2,0]
    hb = ax.hexbin(unit_dict[unit4]['net'],
                  unit_dict[unit4]['rat'], 
                    gridsize=gridsize, cmap='Blues', mincnt=1,
                   )
    ax.set_title('$pp2n$ (5m)', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # pp2n 36
    ax = arr[2,1]
    hb = ax.hexbin(unit_dict[unit5]['net'],
                  unit_dict[unit5]['rat'], 
                    gridsize=gridsize, cmap='Blues', mincnt=1)
    ax.set_title('$pp2n$ (36m)', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    fig.text(0.5, 0.04,
             'Network distance (meters)',ha='center', fontsize=25)
    fig.text(0.07, 0.5,
             'Ratio of Euclidean distance',
             va='center', rotation='vertical', fontsize=25)
    
    if to_fs:
        fs_label = 'fs_'
    else:
        fs_label = ''
    
    plt.savefig('../results/Test_Tract_Leon_FL/plots/%snet_euc_ch2%s.png' % (fs_label, tag), bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    #----------------------------------------------------------------------------------
    
    # Euclidean - Network WEIGHTED
    if to_fs:
        dist_arrays = spgh.load_pickled(mtx_dir, pickle_name='dist_arrays_to_FireStations_Synthetic')
    else:
        dist_arrays = spgh.load_pickled(mtx_dir, pickle_name='dist_arrays')
    
    pp = 'PointPattern'
    census_bases = [#'WeightedParcels',#####################################################
                    'WeightedParcels',
                    'pc2nTracts'+tag,
                    'pc2nBlockGroups'+tag,
                    'pc2nPopulatedBlocks'+tag]
    allocation_bases = ['va2nTracts'+tag,
                        'va2nBlockGroups'+tag,
                        'va2nPopulatedBlocks'+tag,
                        'pp2n5.0Tracts'+tag,
                        'pp2n5.0BlockGroups'+tag,
                        'pp2n5.0PopulatedBlocks'+tag,
                        'pp2n36.0Tracts'+tag,
                        'pp2n36.0BlockGroups'+tag,
                        'pp2n36.0PopulatedBlocks'+tag]
    
    patterns = {name: spgh.load_pickled(cen_dir,
                                        pickle_name='%s_%s' % (pp, name))\
                                        for name in census_bases}
    
    patterns.update({name: spgh.load_pickled(alc_dir,
                                             pickle_name='%s_%s' % (pp, name))\
                                             for name in allocation_bases})
    rep2npts = {}
    for k,p in patterns.items():
        rep2npts[k] = p.df.shape[0]
    
    # dataframe index
    index = list(patterns['pc2nPopulatedBlocks'+tag].segm2pop.keys())
    # instantiate dataframe
    segm2pop_df = pd.DataFrame(index=index)
    # map each {allocation method: segment population} to dataframe
    for name, pattern in patterns.items():
        segm2pop_df[name] = segm2pop_df.index.map(pattern.segm2pop)
    
    
    
    ########################################################
    method1 = 'euclidean'
    method2 = 'network'
    unit1 = 'pc2nPopulatedBlocks'+tag
    unit2 = 'va2nPopulatedBlocks'+tag
    unit3 = 'WeightedParcels'+tag
    unit4 = 'pp2n5.0PopulatedBlocks'+tag
    unit5 = 'pp2n36.0PopulatedBlocks'+tag
    to = '_to_'
    cmtx = '_DataFrame.csv'
    
    units = [unit1, unit2, unit3, unit4, unit5]
    unit_dict = {}
    
    ############################################
    for unit in units:
        seg2pop = segm2pop_df[unit].values
        seg2pop = seg2pop.reshape(1, seg2pop.shape[0])
        
        euc_matrix = dist_arrays['euclidean'] * seg2pop
        euc_vector = euc_matrix[euc_matrix != 0]
        euc_vector = euc_vector / 1000
        
        net_matrix = dist_arrays['network'] * seg2pop
        net_vector = net_matrix[net_matrix != 0]
        net_vector = net_vector / 1000
        
        # calc ratio
        rat_vector = net_vector / euc_vector
        total_size = rat_vector.shape[0]
        
        rat_vector, mask = reject_outliers(rat_vector)
        net_vector = net_vector[mask]
        
        std3 = net_vector.shape[0]
        
        # add to dict
        unit_dict[unit] = {'net': net_vector,
                           'rat': rat_vector,
                           'OD': total_size,
                           'std3': std3}
    
    gridsize = 25
    
    fig, arr = plt.subplots(3,2,figsize=(15,15))
    
    # pc2n
    ax = arr[0,0]
    hb = ax.hexbin(unit_dict[unit1]['net'],
                  unit_dict[unit1]['rat'],
                    gridsize=gridsize,
                   cmap='Greens', mincnt=1)
    ax.set_title('$pc2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # va2n
    ax = arr[0,1]
    hb = ax.hexbin(unit_dict[unit2]['net'],
                  unit_dict[unit2]['rat'], 
                    gridsize=gridsize,
                   cmap='Reds', mincnt=1)
    ax.set_title('$va2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # parcels
    ax = arr[1,0]
    hb = ax.hexbin(unit_dict[unit3]['net'],
                  unit_dict[unit3]['rat'], 
                    gridsize=gridsize, 
                   cmap='RdPu', mincnt=1)
    ax.set_title('$wp2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    if to_fs:
        arr[1,1].text(0, .9, 
                      '$fs$ $\\rightarrow$ $pc2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit1]['std3'], unit_dict[unit1]['OD']), fontsize=18)
        arr[1,1].text(0, .68, 
                      '$fs$ $\\rightarrow$ $va2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit2]['std3'], unit_dict[unit2]['OD']), fontsize=18)
        arr[1,1].text(0, .46, 
                      '$fs$ $\\rightarrow$ $wp2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit3]['std3'], unit_dict[unit3]['OD']), fontsize=18)
        arr[1,1].text(0, .22, 
                      '$fs$ $\\rightarrow$ $pp2n$ (5m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit4]['std3'], unit_dict[unit4]['OD']), fontsize=18)
        arr[1,1].text(0, 0, 
                      '$fs$ $\\rightarrow$ $pp2n$ (36m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit5]['std3'], unit_dict[unit5]['OD']), fontsize=18)
        arr[1,1].axis('off')
    else:
        arr[1,1].text(0, .9, 
                      '$pc2n$ $\\rightarrow$ $pc2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit1]['std3'], unit_dict[unit1]['OD']), fontsize=18)
        arr[1,1].text(0, .68, 
                      '$va2n$ $\\rightarrow$ $va2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit2]['std3'], unit_dict[unit2]['OD']), fontsize=18)
        arr[1,1].text(0, .46, 
                      '$wp2n$ $\\rightarrow$ $wp2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit3]['std3'], unit_dict[unit3]['OD']), fontsize=18)
        arr[1,1].text(0, .22, 
                      '$pp2n$ $\\rightarrow$ $pp2n$ (5m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit4]['std3'], unit_dict[unit4]['OD']), fontsize=18)
        arr[1,1].text(0, 0, 
                      '$pp2n$ $\\rightarrow$ $pp2n$ (36m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit5]['std3'], unit_dict[unit5]['OD']), fontsize=18)
        arr[1,1].axis('off')
    
    # pp2n 5
    ax = arr[2,0]
    hb = ax.hexbin(unit_dict[unit4]['net'],
                  unit_dict[unit4]['rat'], 
                    gridsize=gridsize, 
                   cmap='Blues', mincnt=1)
    ax.set_title('$pp2n$ (5m)', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # pp2n 36
    ax = arr[2,1]
    hb = ax.hexbin(unit_dict[unit5]['net'],
                  unit_dict[unit5]['rat'], 
                    gridsize=gridsize, 
                   cmap='Blues', mincnt=1)
    ax.set_title('$pp2n$ (36m)', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    
    fig.text(0.5, 0.04,
             'Population-weighted network distance (kilometers)',ha='center', fontsize=25)
    fig.text(0.07, 0.5,
             'Ratio of population-weighted Euclidean distance',
             va='center', rotation='vertical', fontsize=25)
    
    if to_fs:
        fs_label = 'fs_'
    else:
        fs_label = ''
    
    plt.savefig('../results/Test_Tract_Leon_FL/plots/%sweighted_net_euc_ch2%s.png' % (fs_label, tag), bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    
    #####################################################################################################
def ch3_plots(adjusted=False, to_fs=True):
    """
    
    """
    
    
    
    
    
    if adjusted:
        tag = '_ADJUSTED'
    else:
        tag = ''
    
    area = 'Leon_FL'
    area_suffix = '_2010.shp'
    place_time_shp = area + area_suffix
    
    data = '%s/%s/clean/' % ('../data', area)
    
    obs_dir = data + 'observation_data/'
    cen_dir = data + 'census_data/'
    net_dir = data + 'network_data/'
    mtx_dir = data + 'cost_matrices/'
    alc_dir = data + 'allocation_data/'
    
    
    #------------------------------------------------------------------------------------------
    # leon tracts and bkgs
    tracts_file = '%sCensusTracts_%s' % (cen_dir, place_time_shp)
    tracts = gpd.read_file(tracts_file)
    bkgs_file = '%sCensusBlockGroups_%s' % (cen_dir, place_time_shp)
    bkgs = gpd.read_file(bkgs_file)
    city_file = '%sTallahassee_%s' % (cen_dir, place_time_shp)
    city = gpd.read_file(city_file)
    
    trtlbl = mpatches.Patch(fc='w', ec='g', alpha=.5, lw=4,
                            label='Tracts\n$\longrightarrow$ $POLY$=%s' % tracts.shape[0])
    bgrlbl = mpatches.Patch(fc='w', ec='b', alpha=.5, lw=2,
                            label='Block Groups\n$\longrightarrow$ $POLY$=%s' % bkgs.shape[0])
    city_handle = mpatches.Patch(fc='gold', ec='gold',
                                      label='Tallahassee city limits',
                                         alpha=.5, linewidth=1)
    
    handles = [trtlbl, bgrlbl, city_handle]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    tracts.plot(ax=ax, color='w',alpha=.5, edgecolor='g', lw=5)
    bkgs.plot(ax=ax, color='w', alpha=.5, edgecolor='b', lw=1)
    city.plot(ax=ax, color='gold', alpha=.5, lw=1, edgecolor='gold', zorder=1)
    
    plt.legend(handles=handles,
                   bbox_to_anchor=(0., -0., 1., -0.),
                   fontsize=13,
                   ncol=3, mode="expand")
    
    x, y = 651000, 144000
    arw = 'rarrow, pad=0.15'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.15, location='lower left')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 578000, 186000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Leon_FL/plots/leon_tracts_bkgs_ch3.png', bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    #------------------------------------------------------------------------------------------
    # leon blocks and unpopulated blocks
    
    blocks_file = '%sCensusBlocks_%s' % (cen_dir, place_time_shp)
    blocks = gpd.read_file(blocks_file)
    
    pop_blocks = blocks[blocks['POP100'] > 0]
    nopop_blocks = blocks[blocks['POP100'] <= 0]
    
    pbklbl = mpatches.Patch(fc='w', ec='k', alpha=.2, ls='-.', lw=1,
                            label='Populated Blocks\n$\longrightarrow$ $POLY$=%s' % pop_blocks.shape[0])
    nbklbl = mpatches.Patch(fc='r', ec='k', alpha=.75, ls='-.',lw=1,
                        label='Unpopulated Blocks\n$\longrightarrow$ $POLY$=%s' % nopop_blocks.shape[0])
    
    handles = [pbklbl, nbklbl]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    pop_blocks.plot(ax=ax, color='w',alpha=.5, edgecolor='k', lw=.25, linestyle='-.')
    nopop_blocks.plot(ax=ax, color='r', alpha=.75, edgecolor='k', lw=.25, linestyle='-.')

    plt.legend(handles=handles,
               bbox_to_anchor=(0., -0, 1., -0),
               fontsize=15,
               ncol=2, mode="expand")
    
    x, y = 651000, 144000
    arw = 'rarrow, pad=0.15'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.15, location='lower left')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 578000, 186000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Leon_FL/plots/leon_blocks_ch3.png', bbox_inches='tight',
                    format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    #------------------------------------------------------------------------------------------
    # leon parcel centroids
    
    centroids_file = '%sWeightedParcels_%s' % (cen_dir,place_time_shp)
    county_centroids = gpd.read_file(centroids_file)
    
    county_no_pop_centroids = county_centroids[county_centroids['SUM_EST_PO'] <= 0.]
    county_pop_centroids = county_centroids[county_centroids['SUM_EST_PO'] > 0.0]
    
    pop_hdl = mlines.Line2D([], [], color='m', markeredgecolor='k',
                            alpha=.5, marker='s', ms=5, linewidth=0,
                                label='Weighted Parcels\n$\\rightarrow$ $RP$=%s' % county_pop_centroids.shape[0]\
                                    + '\n$\\rightarrow$ $\sum pop$=%s' % round(county_pop_centroids['SUM_EST_PO'].sum(), 8))
    
    no_pop_hdl = mlines.Line2D([], [], color='r', markeredgecolor='k',
                            alpha=1., marker='s', ms=6.5, linewidth=0,
                                label='Unweighted Parcels\n$\\rightarrow$ $RP$=%s' % county_no_pop_centroids.shape[0]\
                                       + '\n$\\rightarrow$ $\sum pop$=%s' % county_no_pop_centroids['SUM_EST_PO'].sum())
    
    streets_file = '%sSimplifiedSegms_%s' % (net_dir, place_time_shp)
    streets = gpd.read_file(streets_file)
    street_handle = mlines.Line2D([], [], color='k',
                                  label='Network segments\n$\longrightarrow$ $S$=%s' % streets.shape[0],
                         alpha=.7, linewidth=1)
    
    vertexx_file = '%sSimplifiedNodes_%s' % (net_dir, place_time_shp)
    vertexx = gpd.read_file(vertexx_file)
    vertex_handle = mlines.Line2D([], [], color='k',
                                  label='Network vertices\n$\longrightarrow$ $V$=%s' % vertexx.shape[0],
                                            alpha=1, marker='o', ms=2,
                                            linewidth=0)
    
    handles = [pop_hdl,
           street_handle,
           no_pop_hdl,
           vertex_handle]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    county_pop_centroids.plot(ax=ax, color='m', edgecolor='k', lw=.25,
                       alpha=.5, marker='s', markersize=1.5, zorder=2)
    county_no_pop_centroids.plot(ax=ax, color='r', edgecolor='k', lw=.25,
                       alpha=.75, marker='s', markersize=2., zorder=2)
    streets.plot(ax=ax, color='k', alpha=.5, lw=.35)
    vertexx.plot(ax=ax, color='k', alpha=1, markersize=.35, zorder=1)
    
    plt.legend(handles=handles,
                   bbox_to_anchor=(0., -0., 1., -0.),
                   fontsize=15, borderpad=1., labelspacing=1.,
                   ncol=2, mode="expand", columnspacing=1.)
    
    x, y = 651000, 144000
    arw = 'rarrow, pad=0.15'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.15, location='lower left')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 578000, 186000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Leon_FL/plots/leon_parcel_centroids_ch3.png', bbox_inches='tight',
                    format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    #------------------------------------------------------------------------------------------
    # leon parcels
    
    county_parcel_file = '%sParcelPolygons_%s' % (cen_dir,place_time_shp)
    county_parcels = gpd.read_file(county_parcel_file)
    
    no_tract_id_parcels = county_parcels[county_parcels['TRACT'].isnull()]
    county_no_pop_parcels = county_parcels[(county_parcels['SUM_EST_PO'] <= 0.)\
                     & (~county_parcels['PARCEL_ID'].isin(no_tract_id_parcels['PARCEL_ID']))]
    county_pop_parcels = county_parcels[county_parcels['SUM_EST_PO'] > 0.0]
    
    
    county_ambig_pop_parcels = county_no_pop_parcels[county_no_pop_parcels.PARCEL_ID\
                               .isin(county_no_pop_centroids.PARCEL_ID)]
    
    county_absol_no_pop_parcels = county_no_pop_parcels[~county_no_pop_parcels.PARCEL_ID\
                                      .isin(county_no_pop_centroids.PARCEL_ID)]
    
    pop_hdl = mpatches.Patch(fc='w', ec='k', alpha=1., lw=.1,
                                label='Populated Parcels\n$\\rightarrow$ $POLY$=%s' % county_pop_parcels.shape[0]\
                                    + '\n$\\rightarrow$ $\sum pop$=%s' % round(county_pop_parcels['SUM_EST_PO'].sum(), 8))
    
    ambig_hdl = mpatches.Patch(fc='w', ec='g', alpha=.5, lw=.1, hatch='/////////',
                                label='Ambiguous Parcels\n$\\rightarrow$ $POLY$=%s' % county_ambig_pop_parcels.shape[0]\
                                    + '\n$\\rightarrow$ $\sum pop$=%s' % county_ambig_pop_parcels['SUM_EST_PO'].sum())
    
    no_pop_hdl = mpatches.Patch(fc='r', ec='k', alpha=.75, lw=.1,
                                label='Unpopulated Parcels\n$\\rightarrow$ $POLY$=%s' % county_absol_no_pop_parcels.shape[0]\
                                       + '\n$\\rightarrow$ $\sum pop$=%s' % round(county_absol_no_pop_parcels['SUM_EST_PO'].sum(),20))
    
    no_tract_id_hdl = mpatches.Patch(fc='k', ec='k', alpha=1, lw=2,
                                 label='Excluded Parcels\n$\\rightarrow$ $POLY$=%s' % no_tract_id_parcels.shape[0]\
                                   + '\n$\\rightarrow$ $\sum pop$=%s' % no_tract_id_parcels['SUM_EST_PO'].sum())
    
    handles=[pop_hdl, 
         ambig_hdl,
         no_pop_hdl,
         no_tract_id_hdl]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    county_pop_parcels.plot(ax=ax, color='w', alpha=1, lw=.01, edgecolor='k', zorder=0)
    county_ambig_pop_parcels.plot(ax=ax, color='w', alpha=.5, lw=.01, hatch='////////////',
                        edgecolor='g', zorder=1)
    county_absol_no_pop_parcels.plot(ax=ax, color='r', alpha=.75, lw=.01,
                        edgecolor='k', zorder=1)
    no_tract_id_parcels.plot(ax=ax, color='k', alpha=1,
                             lw=3, edgecolor='k', zorder=3)
    
    plt.legend(handles=handles,
                   bbox_to_anchor=(0., -0., 1., -0.),
                   fontsize=15, borderpad=1., labelspacing=1.,
                   ncol=2, mode="expand", columnspacing=1.)
    
    x, y = 651000, 144000
    arw = 'rarrow, pad=0.15'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.15, location='lower left')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 578000, 186000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Leon_FL/plots/leon_parcels_ch3.png', bbox_inches='tight',
                    format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    #----------------------------------------
    # euclidean - network
    method1 = 'euclidean'
    method2 = 'network'
    unit1 = 'pc2nPopulatedBlocks'
    unit2 = 'va2nPopulatedBlocks'
    unit3 = 'WeightedParcels'
    unit4 = 'pp2n5.0PopulatedBlocks'
    unit5 = 'pp2n36.0PopulatedBlocks'
    to = '_to_'
    cmtx = '_DataFrame.csv'
    
    units = [unit1, unit2, unit3, unit4, unit5]
    unit_dict = {}
    for unit in units:
        
        if to_fs:
            unit_ = 'FireStations'
        else:
            unit_ = unit
        
        # file names
        euc_file = method1 + '_' + unit_ + to + unit + cmtx
        net_file = method2 + '_' + unit_ + to + unit + cmtx
        
        # read files
        net_matrix = pd.read_csv(mtx_dir+net_file, index_col=0)
        net_matrix = net_matrix.values
        
        euc_matrix = pd.read_csv(mtx_dir+euc_file, index_col=0)
        euc_matrix = euc_matrix.values
        
        # calc ratio
        ratio_matrix = net_matrix / euc_matrix
        
        # create vectors
        net_vector = net_matrix.flatten()
        rat_vector = ratio_matrix.flatten()
        
        total_size = rat_vector.shape[0]
        
        # remove outliers more than 1 std deviation away from mean
        rat_vector, mask = reject_outliers(rat_vector)
        net_vector = net_vector[mask]
        
        std3 = net_vector.shape[0]
        
        # add to dict
        unit_dict[unit] = {'net': net_vector,
                            'rat': rat_vector,
                            'OD': total_size,
                            'std3': std3}
                       
    
    gridsize = 25
    
    fig, arr = plt.subplots(3,2,figsize=(15,15),
                            sharex=True, sharey=True)
    
    # pc2n
    ax = arr[0,0]
    hb = ax.hexbin(unit_dict[unit1]['net'],
                  unit_dict[unit1]['rat'],
                    gridsize=gridsize, cmap='Greens', mincnt=1)
    ax.set_title('$pc2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # va2n
    ax = arr[0,1]
    hb = ax.hexbin(unit_dict[unit2]['net'],
                  unit_dict[unit2]['rat'], 
                    gridsize=gridsize, cmap='Reds', mincnt=1)
    ax.set_title('$va2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # parcels
    ax = arr[1,0]
    hb = ax.hexbin(unit_dict[unit3]['net'],
                  unit_dict[unit3]['rat'], 
                    gridsize=gridsize, cmap='RdPu', mincnt=1)
    ax.set_title('$wp2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    
    if to_fs:
        arr[1,1].text(-1, 1.85, 
                      '$fs$ $\\rightarrow$ $pc2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit1]['std3'], unit_dict[unit1]['OD']), fontsize=18)
        arr[1,1].text(-1, 1.63, 
                      '$fs$ $\\rightarrow$ $va2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit2]['std3'], unit_dict[unit2]['OD']), fontsize=18)
        arr[1,1].text(-1, 1.41, 
                      '$fs$ $\\rightarrow$ $wp2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit3]['std3'], unit_dict[unit3]['OD']), fontsize=18)
        arr[1,1].text(-1, 1.19, 
                      '$fs$ $\\rightarrow$ $pp2n$ (5m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit4]['std3'], unit_dict[unit4]['OD']), fontsize=18)
        arr[1,1].text(-1, .97, 
                      '$fs$ $\\rightarrow$ $pp2n$ (36m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit5]['std3'], unit_dict[unit5]['OD']), fontsize=18)
        arr[1,1].axis('off')
    else:
        arr[1,1].text(0, .9, 
                      '$pc2n$ $\\rightarrow$ $pc2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit1]['std3'], unit_dict[unit1]['OD']), fontsize=18)
        arr[1,1].text(0, .68, 
                      '$va2n$ $\\rightarrow$ $va2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit2]['std3'], unit_dict[unit2]['OD']), fontsize=18)
        arr[1,1].text(0, .46, 
                      '$wp2n$ $\\rightarrow$ $wp2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit3]['std3'], unit_dict[unit3]['OD']), fontsize=18)
        arr[1,1].text(0, .22, 
                      '$pp2n$ $\\rightarrow$ $pp2n$ (5m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit4]['std3'], unit_dict[unit4]['OD']), fontsize=18)
        arr[1,1].text(0, 0, 
                      '$pp2n$ $\\rightarrow$ $pp2n$ (36m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit5]['std3'], unit_dict[unit5]['OD']), fontsize=18)
        arr[1,1].axis('off')
    
    # pp2n 5
    ax = arr[2,0]
    hb = ax.hexbin(unit_dict[unit4]['net'],
                  unit_dict[unit4]['rat'], 
                    gridsize=gridsize, cmap='Blues', mincnt=1,
                   )
    ax.set_title('$pp2n$ (5m)', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # pp2n 36
    ax = arr[2,1]
    hb = ax.hexbin(unit_dict[unit5]['net'],
                  unit_dict[unit5]['rat'], 
                    gridsize=gridsize, cmap='Blues', mincnt=1)
    ax.set_title('$pp2n$ (36m)', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    fig.text(0.5, 0.04,
             'Network distance (meters)',ha='center', fontsize=25)
    fig.text(0.07, 0.5,
             'Ratio of Euclidean distance',
             va='center', rotation='vertical', fontsize=25)
    
    if to_fs:
        fs_label = 'fs_'
    else:
        fs_label = ''
    
    plt.savefig('../results/Leon_FL/plots/public_%snet_euc_ch3%s.png' % (fs_label, tag), bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    #----------------------------------------------------------------------------------
    
    # Euclidean - Network WEIGHTED
    if to_fs:
        dist_arrays = spgh.load_pickled(mtx_dir, pickle_name='dist_arrays_to_FireStations')
    else:
        dist_arrays = spgh.load_pickled(mtx_dir, pickle_name='dist_arrays')
    
    pp = 'PointPattern'
    census_bases = [#'WeightedParcels',#####################################################
                    'WeightedParcels',
                    'pc2nTracts'+tag,
                    'pc2nBlockGroups'+tag,
                    'pc2nPopulatedBlocks'+tag]
    allocation_bases = ['va2nTracts'+tag,
                        'va2nBlockGroups'+tag,
                        'va2nPopulatedBlocks'+tag,
                        'pp2n5.0Tracts'+tag,
                        'pp2n5.0BlockGroups'+tag,
                        'pp2n5.0PopulatedBlocks'+tag,
                        'pp2n36.0Tracts'+tag,
                        'pp2n36.0BlockGroups'+tag,
                        'pp2n36.0PopulatedBlocks'+tag]
    
    patterns = {name: spgh.load_pickled(cen_dir,
                                        pickle_name='%s_%s' % (pp, name))\
                                        for name in census_bases}
    
    patterns.update({name: spgh.load_pickled(alc_dir,
                                             pickle_name='%s_%s' % (pp, name))\
                                             for name in allocation_bases})
    rep2npts = {}
    for k,p in patterns.items():
        rep2npts[k] = p.df.shape[0]
    
    # dataframe index
    index = list(patterns['pc2nPopulatedBlocks'+tag].segm2pop.keys())
    # instantiate dataframe
    segm2pop_df = pd.DataFrame(index=index)
    # map each {allocation method: segment population} to dataframe
    for name, pattern in patterns.items():
        segm2pop_df[name] = segm2pop_df.index.map(pattern.segm2pop)
    
    ########################################################
    method1 = 'euclidean'
    method2 = 'network'
    unit1 = 'pc2nPopulatedBlocks'+tag
    unit2 = 'va2nPopulatedBlocks'+tag
    unit3 = 'WeightedParcels'+tag
    unit4 = 'pp2n5.0PopulatedBlocks'+tag
    unit5 = 'pp2n36.0PopulatedBlocks'+tag
    to = '_to_'
    cmtx = '_DataFrame.csv'
    
    units = [unit1, unit2, unit3, unit4, unit5]
    unit_dict = {}
    
    ############################################
    for unit in units:
        seg2pop = segm2pop_df[unit].values
        #seg2pop = seg2pop.reshape(seg2pop.shape[0], 1)
        seg2pop = seg2pop.reshape(1, seg2pop.shape[0])
        
        euc_matrix = dist_arrays['euclidean'] * seg2pop
        euc_vector = euc_matrix[euc_matrix != 0]
        euc_vector = euc_vector / 1000
        
        net_matrix = dist_arrays['network'] * seg2pop
        net_vector = net_matrix[net_matrix != 0]
        net_vector = net_vector / 1000
        
        # calc ratio
        rat_vector = net_vector / euc_vector
        total_size = rat_vector.shape[0]
        
        rat_vector, mask = reject_outliers(rat_vector)
        net_vector = net_vector[mask]
        
        std3 = net_vector.shape[0]
        
        # add to dict
        unit_dict[unit] = {'net': net_vector,
                           'rat': rat_vector,
                           'OD': total_size,
                           'std3': std3}
    
    gridsize = 25
    
    fig, arr = plt.subplots(3,2,figsize=(15,15))
    
    # pc2n
    ax = arr[0,0]
    hb = ax.hexbin(unit_dict[unit1]['net'],
                  unit_dict[unit1]['rat'],
                    gridsize=gridsize,
                   cmap='Greens', mincnt=1)
    ax.set_title('$pc2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # va2n
    ax = arr[0,1]
    hb = ax.hexbin(unit_dict[unit2]['net'],
                  unit_dict[unit2]['rat'], 
                    gridsize=gridsize,
                   cmap='Reds', mincnt=1)
    ax.set_title('$va2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # parcels
    ax = arr[1,0]
    hb = ax.hexbin(unit_dict[unit3]['net'],
                  unit_dict[unit3]['rat'], 
                    gridsize=gridsize, 
                   cmap='RdPu', mincnt=1)
    ax.set_title('$wp2n$', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    if to_fs:
        arr[1,1].text(0, .9, 
                      '$fs$ $\\rightarrow$ $pc2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit1]['std3'], unit_dict[unit1]['OD']), fontsize=18)
        arr[1,1].text(0, .68, 
                      '$fs$ $\\rightarrow$ $va2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit2]['std3'], unit_dict[unit2]['OD']), fontsize=18)
        arr[1,1].text(0, .46, 
                      '$fs$ $\\rightarrow$ $wp2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit3]['std3'], unit_dict[unit3]['OD']), fontsize=18)
        arr[1,1].text(0, .22, 
                      '$fs$ $\\rightarrow$ $pp2n$ (5m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit4]['std3'], unit_dict[unit4]['OD']), fontsize=18)
        arr[1,1].text(0, 0, 
                      '$fs$ $\\rightarrow$ $pp2n$ (36m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit5]['std3'], unit_dict[unit5]['OD']), fontsize=18)
        arr[1,1].axis('off')
    else:
        arr[1,1].text(0, .9, 
                      '$pc2n$ $\\rightarrow$ $pc2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit1]['std3'], unit_dict[unit1]['OD']), fontsize=18)
        arr[1,1].text(0, .68, 
                      '$va2n$ $\\rightarrow$ $va2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit2]['std3'], unit_dict[unit2]['OD']), fontsize=18)
        arr[1,1].text(0, .46, 
                      '$wp2n$ $\\rightarrow$ $wp2n$\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit3]['std3'], unit_dict[unit3]['OD']), fontsize=18)
        arr[1,1].text(0, .22, 
                      '$pp2n$ $\\rightarrow$ $pp2n$ (5m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit4]['std3'], unit_dict[unit4]['OD']), fontsize=18)
        arr[1,1].text(0, 0, 
                      '$pp2n$ $\\rightarrow$ $pp2n$ (36m)\n$\quad$ $OD^\prime$=%s ($OD$=%s)'\
                      % (unit_dict[unit5]['std3'], unit_dict[unit5]['OD']), fontsize=18)
        arr[1,1].axis('off')
    
    # pp2n 5
    ax = arr[2,0]
    hb = ax.hexbin(unit_dict[unit4]['net'],
                  unit_dict[unit4]['rat'], 
                    gridsize=gridsize, 
                   cmap='Blues', mincnt=1)
    ax.set_title('$pp2n$ (5m)', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    # pp2n 36
    ax = arr[2,1]
    hb = ax.hexbin(unit_dict[unit5]['net'],
                  unit_dict[unit5]['rat'], 
                    gridsize=gridsize, 
                   cmap='Blues', mincnt=1)
    ax.set_title('$pp2n$ (36m)', fontsize=20)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('OD density')
    
    
    fig.text(0.5, 0.04,
             'Population-weighted network distance (kilometers)',ha='center', fontsize=25)
    fig.text(0.07, 0.5,
             'Ratio of population-weighted Euclidean distance',
             va='center', rotation='vertical', fontsize=25)
    
    if to_fs:
        fs_label = 'fs_'
    else:
        fs_label = ''
    
    plt.savefig('../results/Leon_FL/plots/public_%sweighted_net_euc_ch3%s.png' % (fs_label, tag), bbox_inches='tight',
                format='png', dpi=400 ,quality=95)
    
    plt.close()
    
    
    
def ch4_plots(adjusted=False, to_fs=True):
    """
    
    """
    
    if adjusted:
        tag = '_ADJUSTED'
    else:
        tag = ''
    
    area = 'Leon_FL'
    area_suffix = '_2010.shp'
    place_time_shp = area + area_suffix
    
    
    clean_data = '../%s/%s/clean/' % ('data', area)
    init_data = '../%s/%s/initial/' % ('data', area)
    
    clean_obs_dir = clean_data + 'observation_data/'
    init_obs_dir = init_data + 'FireStations_Leon_FL/'
    
    cen_dir = clean_data + 'census_data/'
    net_dir = clean_data + 'network_data/'
    
    ##########
    clean_fire_stations_file = clean_obs_dir + 'FireStations_Leon_FL_2010.shp'
    clean_fire_stations = gpd.read_file(clean_fire_stations_file)
    
    ##########
    initial_fire_stations_file = init_obs_dir + 'FireStations_Leon_FL.shp'
    initial_fire_stations = gpd.read_file(initial_fire_stations_file)
    #out_subset_cols = ['NAME', 'ADDRESS']
    #out_subset = initial_fire_stations[out_subset_cols]
    #out_subset['Permanent_2010'] = ['Yes'] * 15 + ['No'] * 6
    
    #######
    initial_fire_stations = initial_fire_stations.to_crs(clean_fire_stations.crs)
    
    ##
    initial_only = set(initial_fire_stations.ADDRESS)\
                   .symmetric_difference(\
                   set(clean_fire_stations.ADDRESS))
    
    ###
    initial_only = initial_fire_stations[\
                   initial_fire_stations.ADDRESS.isin(initial_only)]
    
    ###
    county_file = '%sCounties_%s' % (cen_dir, place_time_shp)
    county = gpd.read_file(county_file)
    
    ###
    
    streets_file = '%sSimplifiedSegms_%s' % (net_dir, place_time_shp)
    streets = gpd.read_file(streets_file)
    
    ###
    fs_clean_handle = mlines.Line2D([], [], color='r', marker='$F$',
                            label='Fire stations included\n$\longrightarrow$ $F$=%s' % clean_fire_stations.shape[0],
                            alpha=1, ms=15, linewidth=0)
    
    fs_initial_handle = mlines.Line2D([], [], color='k', marker='$f$',
                                label='Fire stations excluded\n$\longrightarrow$ $f$=%s' % initial_only.shape[0],
                                alpha=1, ms=15, linewidth=0)
    
    street_handle = mlines.Line2D([], [], color='k',
                                  label='Network segments',
                                  alpha=.25, linewidth=.5)
    
    county_handle = mlines.Line2D([], [], color='k',
                                  label='Leon County boundary',
                                     alpha=.75, linewidth=1)
    
    handles = [fs_clean_handle, street_handle, fs_initial_handle, county_handle]
    
    fig, ax = plt.subplots(figsize=(9, 9))
    county.plot(ax=ax, color='w', alpha=.75, lw=1, edgecolor='k', zorder=0)
    streets.plot(ax=ax, color='k', alpha=.25, lw=.5, zorder=0)
    initial_only.plot(ax=ax, color='b', edgecolor='k',
                           alpha=1, marker='$f$', markersize=50, zorder=2)
    clean_fire_stations.plot(ax=ax, color='r', alpha=1, marker='$F$', 
                             markersize=50, zorder=2)
    
    plt.legend(handles=handles,
               bbox_to_anchor=(0., -0, 1., -0),
               fontsize=15,
               ncol=2, mode="expand")
    
    x, y = 651000, 144000
    arw = 'rarrow, pad=0.15'
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    ax.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)
    
    scalebar = ScaleBar(1, units='m', frameon=False, pad=.5,
                        font_properties={'size':8, 'style':'italic'},
                        length_fraction=.15, location='lower left')
    ax.add_artist(scalebar)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    proj = 'epsg 2779: NAD83(HARN)/Florida North'
    x, y = 578000, 186000
    ax.text(x, y, proj, fontsize=6.5,
              fontstyle='italic', bbox=bbox_props);
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('../results/Leon_FL/plots/leon_fs_ch4.png', bbox_inches='tight',
                    format='png', dpi=400 ,quality=95)
   
    plt.close()
    
    
    