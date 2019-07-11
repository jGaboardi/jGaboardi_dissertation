"""hard coded plotting of dissertation figures
"""

# standard library imports
import sys

# non-standard library imports
import geopandas as gpd
import pandas as pd
import numpy as np


# project library imports
from at1866_Master import spaghetti as spgh


def ch2_table1():
    """
    pp2n : {float, int}
        pp2n parallel offset parameter
    """
    
    area_prefix = 'Test_Grid_'
    area_base = 'Leon_FL'
    area_suffix = '_2010.shp'
    area = area_prefix + area_base
    place_time_shp = area_base + area_suffix
    
    data = '%s/%s/clean/' % ('../data', area)
    
    obs_dir = data + 'observation_data/'
    cen_dir = data + 'census_data/'
    net_dir = data + 'network_data/'
    mtx_dir = data + 'cost_matrices/'
    alc_dir = data + 'allocation_data/'
    
    pp = 'PointPattern'
    census_bases = ['HouseholdsSynthetic',
                    'pc2nPopulatedBlocks',]
    allocation_bases = ['va2nPopulatedBlocks',
                        'pp2n0.25PopulatedBlocks']

    patterns = {name: spgh.load_pickled(cen_dir,
                                        pickle_name='%s_%s' % (pp, name))\
                                        for name in census_bases}

    patterns.update({name: spgh.load_pickled(alc_dir,
                                             pickle_name='%s_%s' % (pp, name))\
                                             for name in allocation_bases})

    # dataframe index
    index = list(patterns['pc2nPopulatedBlocks'].segm2pop.keys())

    # instantiate dataframe
    segm2pop_df = pd.DataFrame(index=index)

    # map each {allocation method: segment population} to dataframe
    for name, pattern in patterns.items():
        segm2pop_df[name] = segm2pop_df.index.map(pattern.segm2pop)
        if not name.startswith('H'):
            T = segm2pop_df['HouseholdsSynthetic']
            E = segm2pop_df[name]
            col = '$e(%s)$' % name[:4]
            segm2pop_df[col] = T - E

    # name index
    segm2pop_df.index.name = patterns['pc2nPopulatedBlocks'].sid_name

    sums = {c:segm2pop_df[c].sum() for c in segm2pop_df.columns}

    segm2pop_df = segm2pop_df.append(sums,ignore_index=True)

    segm2pop_df.index.name = 'SegID'

    full_names = ['HouseholdsSynthetic',
                  'pc2nPopulatedBlocks',
                  'va2nPopulatedBlocks',
                  'pp2n0.25PopulatedBlocks',]
    new_names = ['Synth. Households','$pc2n$', '$va2n$', '$pp2n$',]
    rename_cols = dict(zip(full_names, new_names))

    segm2pop_df.rename(index=str, columns=rename_cols, inplace=True)

    segm2pop_df.to_csv('../results/Test_Grid_Leon_FL/tables/ch2_table1.csv',
                       float_format='%.6f')



def ch2_table2_table3(area, adjusted=False, chapter=None):
    """
    """
    
    if adjusted:
        tag = '_ADJUSTED'
    else:
        tag = ''
    
    
    #area_suffix = '_2010' + tag
    area_suffix = tag
    area = area + area_suffix
    
    df1_out = '../results/%s/tables/ch%s_table2%s.csv' % (area, chapter, tag)
    df2_out = '../results/%s/tables/ch%s_table3%s.csv' % (area, chapter, tag)
    
    data = '%s/%s/clean/' % ('../data', area)
    
    obs_dir = data + 'observation_data/'
    cen_dir = data + 'census_data/'
    net_dir = data + 'network_data/'
    mtx_dir = data + 'cost_matrices/'
    alc_dir = data + 'allocation_data/'

    pp = 'PointPattern'
    census_bases = [#'WeightedParcels',############################################### HOUSHOLDS
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
    # name index
    segm2pop_df.index.name = patterns['pc2nPopulatedBlocks'+tag].sid_name

    non_stat_cols = segm2pop_df.columns
    segm2pop_df['Mean'] = segm2pop_df[non_stat_cols].mean(axis=1)
    segm2pop_df['StD'] = segm2pop_df[non_stat_cols].std(axis=1)

    descriptive_stats = ['RP', 'sum', 'perc_0_segms', 'min', 'min_wo_0', 'max',
                         'mean', 'mean_wo_0', 'std', 'std_wo_0']

    desc_stat_df = pd.DataFrame(columns=descriptive_stats,
                                index=segm2pop_df.columns[:-2])

    for idx in desc_stat_df.index:
        RP = rep2npts[idx]
        col = segm2pop_df[idx]
        tot_segms = col.shape[0]
        zero_segms = tot_segms - col[col > 0].shape[0]
        perc_zero = (float(zero_segms) / float(tot_segms)) * 100.00

        stats = [RP,
                 col.sum(),
                 perc_zero,
                 col.min(),
                 col[col > 0].min(),
                 col.max(),
                 col.mean(),
                 col[col > 0].mean(),
                 col.std(),
                 col[col > 0].std()]

        desc_stat_df.loc[idx] = stats

    desc_stat_df = desc_stat_df.astype(float)

    new_names = ['wp2n', 'pc2nTracts', 'pc2nBlockGroups', 'pc2nPopulatedBlocks', 
            'va2nTracts', 'va2nBlockGroups', 'va2nPopulatedBlocks', 
            'pp2n5.0Tracts', 'pp2n5.0BlockGroups', 'pp2n5.0PopulatedBlocks', 
            'pp2n36.0Tracts', 'pp2n36.0BlockGroups', 'pp2n36.0PopulatedBlocks']

    desc_stat_df.index = new_names

    new_keys = ['pc2nTracts', 'va2nTracts', 'pp2n5.0Tracts','pp2n36.0Tracts',
            'pc2nBlockGroups', 'va2nBlockGroups', 'pp2n5.0BlockGroups',
            'pp2n36.0BlockGroups', 'pc2nPopulatedBlocks', 'va2nPopulatedBlocks',
            'pp2n5.0PopulatedBlocks', 'pp2n36.0PopulatedBlocks', 'wp2n']

    new_df = pd.DataFrame(columns=desc_stat_df.columns, index=new_keys)

    for idx in new_df.index:
        new_df.loc[idx] = desc_stat_df.loc[idx]
    
    df1 = new_df[['RP', 'perc_0_segms', 'sum', 'max']]
    df1.index.name = 'Method'
    df1.to_csv(df1_out, float_format='%.6f')
                            
    df2 = new_df[['min', 'mean', 'std', 'min_wo_0','mean_wo_0', 'std_wo_0']]
    df1.index.name = 'Method'
    df2.to_csv(df2_out, float_format='%.6f')
