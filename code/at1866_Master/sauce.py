"""
functions and utilities for spaghetti
this include functions used for cleaning both tiger roads
and tiger edges files
"""

# standard library imports
import os, re, zipfile, copy, io
try:
    import requests
except:
    print('Could not import package: requests')
from ast import literal_eval

# non-standard library imports
import geopandas as gpd
import numpy as np
from shapely.geometry import Point, MultiPoint
from shapely.geometry import LineString, MultiLineString
from shapely.geometry import GeometryCollection
from shapely.ops import linemerge, polygonize

# used to supress warning in addIDX()
gpd.pd.set_option('mode.chained_assignment',None)

# project library imports
from . import utils


###############################################################################
################ TIGER/Line clean up functionality ############################
###############################################################################

def get_raw_tiger_edges(net, tiger_type='Edges'):
    """set paths for raw tiger data
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    tiger_type : str
        current fucntionality supports 'Edges'
    
    Returns
    -------
    shp_file : str
        update 'initial' path to point to .shp
    """
    
    # descriptive file name
    shp_desc = tiger_type + net.place_time
    
    # upper directory to .shp
    shp_path = net.segmdata + shp_desc
    
    # full path to local .shp including ,shp itself
    shp_file = shp_path + '/' + shp_desc + net.file_type
    
    # if the file already exists return the file path
    if os.path.exists(shp_file):
        return shp_file
    co_st = net.county + '_' + net.state
    
    # for test areas
    if net.study_area != co_st:
        
        # path to full set of data
        full_file = net.segmdata.replace(net.study_area, co_st)
        
        # make sure full county exits
        full_shp_file = full_file + shp_desc + '/' + shp_desc + net.file_type
        if not os.path.exists(full_shp_file):
            raise Exception(full_shp_file + ' does not exist.')
        
        # make subset directory if it doesn't exist
        initial_dir = net.segmdata + shp_desc
        if not os.path.exists(initial_dir):
            os.makedirs(initial_dir)
        
        # subset roads from full file
        subset_roads = subset_for_tests(net.study_area, net.proj_trans)
        subset_roads.to_file(shp_file)
        
        return shp_file
    
    # for full counties
    else:
        
        # FIPS dictionary
        st_fp, ct_fp = utils.get_fips(net.state, net.county)
        
        # file name on site
        remote_name = 'tl_' + net.year + '_'\
                      + st_fp + ct_fp + '_' + tiger_type.lower()
        remote_dir = 'geo/tiger/TIGER' + net.year\
                     + '/' + tiger_type.upper() + '/'
        
        # zip file name to be written
        zip_file_name = shp_path + '.zip'
        
        # stream from url -- set url for zipped file location
        zip_file_url = 'https://www2.census.gov/'\
                       + remote_dir + remote_name + '.zip'
        
        # send request
        r = requests.get(zip_file_url, stream=True)
        
        # isolate the zipped file and extract
        zip_file = zipfile.ZipFile(io.BytesIO(r.content))
        
        # write file to disk
        zip_file.extractall(shp_path)
        zip_file.close()
        
        # rename files descriptively
        for filename in os.listdir(shp_path):
            
            # split on '.' to handle .shp.xml
            filename_split = re.split(r'[\.]',filename)
            filename_split[0] = shp_desc
            
            # descriptive name
            desc_name = '.'.join(filename_split)
            rename_from = shp_path + '/' + filename
            rename_to = shp_path + '/' + desc_name
            
            # rename files
            os.rename(rename_from, rename_to)
        
        return shp_file


def subset_for_tests(subset, crs):
    """Extract network segments intersecting the tract buffer
    
    Parameters
    ----------
    subset : str
        network subset name
    crs : int
        projected coordinate reference system
    
    Returns
    -------
    tract_segms : geopandas.GeoDataFrame
        network segments intersecting the tract buffer
    """
    
    if subset == 'Test_Grid_Leon_County':
        raise Exception('does not exist')
    
    else:
        
        # read in roads
        full_roads_file =\
        '../data/Leon_FL/clean/network_data/SimplifiedSegms_Leon_FL_2010.shp'
        
        full_roads = gpd.read_file(full_roads_file)
        full_roads = utils.set_crs(full_roads, proj_init=crs)
        
        # read in tracts
        hard_code = '../data/' + subset\
                        + '/clean/census_data/CensusTracts_Leon_FL_2010.shp'
        tract = gpd.read_file(hard_code)
        tract = utils.set_crs(tract, proj_init=crs)
        
        # extract geometry and create small buffer to get border roads
        tract_ = tract.geometry.squeeze().buffer(1)
        
        # network segments intersecting the tract buffer
        tract_segms = full_roads[full_roads.intersects(tract_)]
        tract_segms = utils.set_crs(tract_segms, proj_init=crs)
        
        return tract_segms


def tiger_netprep(net, in_file=None, calc_len=False):
    """scrub a raw tiger lines file at the county level to prep for
    network. either TIGER/Line EDGES or ROADS (preferably EDGES)
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    in_file : str
        file path to raw data. Default is None.
    calc_len : bool
        calculated length and add column. Default is False.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        fully cleansed roads ready to process network instantiation
    """
    
    # Reproject roads and subset by road type
    subphase_file = net.inter+'StreetsSubset_Initial'
    full_file = net.inter+'StreetsFull'
    gdf = initial_subset(net, in_file, calc_len=calc_len,
                         save_full=full_file, save_subphase=subphase_file)
    
    # Correcting ring roads
    subphase_file = net.inter+'RingRoadsCorrected'
    gdf = ring_correction(net, gdf, save_keep=subphase_file)
    
    # Cleanse SuperCycle
    # Before splitting:
    # TIGER EDGES - Weld interstate segment pieces
    # TIGER ROADS - Drop equal & contained geoms in sequence
    gdf = cleanse_supercycle(net, gdf, series=True, inherit_attrs=True,
                             calc_len=calc_len, equal_geom=True,
                             contained_geom=True)
    
    return gdf


def initial_subset(net, raw_data, calc_len=False,
                   save_full=None, save_subphase=None):
    """Initalize a network data cleanse from raw tiger line files
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    raw_data : str
        directory to find raw tiger data
    calc_len : bool
        calculated length and add column. Default is False.
    save_full : str
        directory to save full files. Default is None.
    save_subphase : str
        path to save subphase file Default is None.
    
    Returns
    -------
    gdf  : geopandas.GeoDataFrame
        initial dataframe of scrubbed roads
    """
    
    # Read in raw TIGER street data
    gdf = setup_raw(net, raw_data=raw_data, calc_len=calc_len)
    
    # remove trouble maker segments
    if net.discard_segs:
        gdf = gdf[~gdf[net.attr2].isin(net.discard_segs)]
    
    # Add three new MTFCC columns for feature class,
    # description, and rank
    if net.mtfcc_types:
        mtfcc_columns = ['FClass', 'Desc', 'MTFCCRank']
        for c in mtfcc_columns:
            gdf[c] = [net.mtfcc_types[mtfcc_type][c]\
                      for mtfcc_type in gdf[net.attr1]]
    
    if save_full:
        gdf.to_file(save_full+net.file_type)
    
    # Subset roads
    if 'FClass' in gdf.columns and net.mtfcc_discard:
        gdf = utils.record_filter(gdf, column='FClass',
                                  mval=net.mtfcc_discard, oper='out')
    gdf.reset_index(drop=True, inplace=True)
    gdf = label_rings(gdf, geo_col=net.geo_col)
    
    # create segment xyID
    segm2xyid = utils.generate_xyid(df=gdf, geom_type='segm',
                                    geo_col=net.geo_col)
    gdf = utils.fill_frame(gdf, col=net.xyid, data=segm2xyid)
    
    if save_subphase:
        gdf.to_file(save_subphase+net.file_type)
    
    return gdf


def label_rings(df, geo_col=None):
    """label each line segment as ring (True) or not (False)
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        dataframe of geometries
    geo_col : str
        geometry column name. Default is None.
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        dataframe of geometries
    """
    
    df['ring'] = ['False']*df.shape[0]
    for idx in df.index:
        if df[geo_col][idx].is_ring:
            df['ring'][idx] = 'True'
    
    return df


def setup_raw(net, raw_data=None, calc_len=False):
    """initial raw data prep
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    raw_data : str
        directory to find raw tiger data
    calc_len : bool
        calculated length and add column. Default is False.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        initial dataframe of scrubed roads
    """
    
    gdf = gpd.read_file(raw_data)
    
    try:
        # tiger edges subsets
        if net.tiger_edges and net.edge_subsets:
            for ss in net.edge_subsets:
                gdf = gdf[ss[1]['relate'](gdf[ss[1]['col']], ss[1]['val'])]
    except KeyError:
        pass
    
    # weld multilinestrings segs into linestrings -- only occurs in 'roads'
    if net.tiger_roads:
        gdf = seg_welder(gdf, geo_col=net.geo_col)
    init_records = gdf.shape[0]
    gdf['Phase'] = 'Initial'
    
    if net.proj_trans: # Transform
        gdf = utils.set_crs(gdf, proj_init=net.proj_init,
                            proj_trans=net.proj_trans)
    
    if calc_len:
        gdf = add_length(gdf, len_col=net.len_col, geo_col=net.geo_col)
        initial_len = gdf[net.len_col].sum()
    
    # segment count and length for the raw tiger data
    raw_data_info = {'segment_count':init_records, 'length':initial_len}
    net.raw_data_info = raw_data_info
    
    return gdf


def add_length(frame, len_col=None, geo_col=None):
    """add length column to a dataframe

    Parameters
    ----------
    frame : geopandas.GeoDataFrame
        dataframe of geometries
    len_col : str
        length column name in dataframe. Default is None.
    geo_col : str
        geometry column name. Default is None.
    
    Returns
    -------
    frame : geopandas.GeoDataFrame
        updated dataframe of geometries
    """
    
    if list(frame.columns).__contains__(len_col):
        frame = frame.drop(len_col, axis=1)
    frame[len_col] = [frame[geo_col][idx].length for idx in frame.index]
    
    return frame


def ring_correction(net, df, save_keep=None):
    """ring roads should start and end with the point at which it
    intersects with another road segment. This algorithm find instances
    where rings roads are digitized incorrectly which results in ring
    roads having their endpoints somewhere in the middle of the line,
    then  corrects the loop by updating the geometry. Length and
    attributes of the original line segment are not changed.
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    df : geopandas.GeoDataFrame
        dataframe of road segments
    save_keep  : str
        path to save the shapefile. Default is None.
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        updated dataframe of road segments
    """
    
    # subset only ring roads
    rings_df = df[df['ring'] == 'True']
    ringsidx, corrected_rings = rings_df.index, 0
    
    for idx in ringsidx:
        LOI = rings_df[net.geo_col][idx]
        
        # get idividual ring road - normal road pairs intersection
        i_geoms = get_intersecting_geoms(net, df1=df, geom1=idx, wbool=False)
        i_geoms = i_geoms[i_geoms.index != idx]
        
        # rings that are not connected to the network will be removed
        if i_geoms.shape[0] < 1:
            continue
        node = i_geoms[net.geo_col][:1].intersection(LOI).values[0]
        
        # if pre cleaned and segments still overlap
        if type(node) != Point:
            continue
        
        node_coords = list(zip(node.xy[0], node.xy[1]))
        line_coords = list(zip(LOI.coords.xy[0], LOI.coords.xy[1]))
        
        # if problem ring road
        # (e.g. the endpoint is not the intersection)
        if node_coords[0] != line_coords[0]:
            updated_line = _correct_ring(node_coords, line_coords)
            
            # update dataframe record
            df[net.geo_col][idx] = updated_line
            corrected_rings += 1
    
    df.reset_index(drop=True, inplace=True)
    df = add_ids(df, id_name=net.sid_name)
    
    # add updated xyid
    segm2xyid = utils.generate_xyid(df=df, geom_type='segm',
                                    geo_col=net.geo_col)
    df = utils.fill_frame(df, col=net.xyid, data=segm2xyid)
    
    # corrected ring road count
    net.corrected_rings = corrected_rings
    if save_keep:
        df.to_file(save_keep+net.file_type)
    
    return df


def _correct_ring(node_coords, line_coords, post_check=False):
    """helper function for ring_correction
    
    Parameters
    ----------
    node_coords : list
        xy tuple for a node.
    line_coords  : list
        all xy tuple for a line.
    post_check : bool
        check following a cleanse cycle. Default is False.
    
    Returns
    -------
    updated_line : shapely.LineString
        ring road updated so that it begins and ends at
        the intersecting node.
    """
    
    if post_check:
        node_coords = [node_coords]
        if node_coords[0] == line_coords[0]:
            return line_coords
    
    for itemidx, coord in enumerate(line_coords):
        # find the index of the intersecting coord in the line
        if coord == node_coords[0]:
            break
    
    # adjust the line coordinates for the true ring start/end
    updated_line = line_coords[itemidx:] + line_coords[1:itemidx+1]
    updated_line = LineString(updated_line)
    
    return updated_line


def get_intersecting_geoms(net, df1=None, geom1=None,
                           df2=None, geom2=None, wbool=True):
    """return the subset of intersecting geometries
    from within a geodataframe
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    df1 : geopandas.GeoDataFrame
        primary dataframe. Default is None.
    geom1 : int
        geometry index. Default is None.
    df2 : geopandas.GeoDataFrame
        secondary dataframe . Default is None.
    geom2 : int
        geometry index. Default is None.
    wbool : bool
        return a boolean object for intersections. Default is True.
    geo_col : str
        geometry column name. Default is None.
    
    Returns
    -------
    i_geom : geopandas.GeoDataFrame
        intersecting geometry subset
    i_bool : bool array
        optional return of intersecting geoms
    """ 
    
    # if there *IS NO* dataframe 2 in play
    if not hasattr(df2, net.geo_col):
        i_bool = df1.intersects(df1[net.geo_col][geom1])
    
    # if there *IS* dataframe 2 in play
    else:
        i_bool = df1.intersects(df2[net.geo_col][geom2])
    i_geom = df1[i_bool]
    
    if wbool:
        return i_bool, i_geom
    
    else:
        return i_geom


def add_ids(frame, id_name=None):
    """add an idx column to a dataframe
    
    Parameters
    ----------
    frame : geopandas.GeoDataFrame
        dataframe of geometries
    id_name : str
        name of id column. Default is None.
    
    Returns
    -------
    frame : geopandas.GeoDataFrame
        updated dataframe of geometries
    """
    
    frame[id_name] = [idx for idx in range(frame.shape[0])]
    frame[id_name] = frame[id_name].astype(int)
    
    return frame


def cleanse_supercycle(net, gdf, series=False, inherit_attrs=False,
                       calc_len=True, equal_geom=False, contained_geom=False,
                       save_out=True, skip_restr=True):
    """One iteration of a cleanse supercycle; then repeat as
    necessary --> 1. Drop equal geoms; 2. Drop contained geoms;
    3. Split line segments
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    gdf : geopandas.GeoDataFrame
        streets dataframe
    inherit_attrs : bool
        inherit attributes from the dominant line segment.
        Default is False.
    series : bool
        search a geoseries. Default is False.
    equal_geom : bool
        drop all but one of equal geometries. Default is False.
    contained_geom : bool
        drop all contained geometries. Default is False.
    calc_len : bool
        calculated length and add column. Default is True.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        streets dataframe
    """
    
    iteration = 0
    completely_scrubbed = False
    while not completely_scrubbed:
        iteration += 1
        itr = str(iteration)
        
        if net.tiger_edges:
            gdf = restriction_welder(net, gdf)
        
        if net.tiger_roads:
            
            # Drop equal geoms
            if net.inter:
                inter_file = net.inter+'NoEqualGeom'+itr
            else:
                inter_file = None
            gdf = streets_filter(net, gdf, series=series,
                                 equal_geom=equal_geom, save_keep=inter_file)
            
            # Drop contained geometries
            if net.inter:
                inter_file = net.inter+'NoContainedGeom'+itr
            else:
                inter_file = None
            gdf = streets_filter(net, gdf, series=series, save_keep=inter_file,
                                 contained_geom=contained_geom)
        
        # Split segments
        if net.inter:
            inter_file = net.inter+'SplitSegms'+itr
        else:
            inter_file = None
        
        gdf = line_splitter(net, gdf, proj_init=net.proj_trans,
                            calc_len=calc_len, save_keep=inter_file, stage=itr,
                            inherit_attrs=inherit_attrs)
        
        # Determine if completely scrubbed
        if not net.tiger_edges:
            completely_scrubbed = check_cleanliness(net, gdf)
        else:
            completely_scrubbed = True
    
    # cycles needed for cleanse
    net.cleanse_cycles = iteration
    net.scrubbed = completely_scrubbed
    
    # Re-lablel Rings
    gdf = label_rings(gdf, geo_col=net.geo_col)
    gdf.reset_index(inplace=True, drop=True)
    
    return gdf


def restriction_welder(net, gdf, phase='Phase', phase_name='InterstateWeld'):
    """weld each set of restricted segments (e.g. interstates)
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    gdf : geopandas.GeoDataFrame
        streets dataframe
    phase : str
        column name in dataframe. Default is 'Phase'.
    phase_name : str
        name of phase. Default is 'InterstateWeld'.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        updated streets dataframe
    """
    
    # make restricted subset
    restr_ss = gdf[gdf[net.attr1] == net.mtfcc_split]
    
    try:
        restr_names = [str(grp) for grp\
                                in restr_ss[net.mtfcc_split_grp].unique()]
    except KeyError:
        return gdf
    
    # create a sub-subset for each group (e.g. interstate)
    for grp in restr_names:
        ss = restr_ss[restr_ss[net.mtfcc_split_grp]== grp]
        
        # get restriction segments to restriction nodes lookup dict
        # and restriction nodes to restriction segments lookup dict
        s2n, n2s = associate(initial_weld=True, net=net, df=restr_ss, ss=ss)
        
        # x2x topologies
        s2s = get_neighbors(s2n, n2s, astype=dict)
        n2n = get_neighbors(n2s, s2n, astype=dict)
        
        # get rooted connected components
        s2s_cc = get_roots(s2s)
        
        # weld together segments from each component of the group
        for cc in s2s_cc:
            keep_id, all_ids = cc[0], cc[1]
            drop_ids = copy.deepcopy(all_ids)
            drop_ids.remove(keep_id)
            
            # subset of specific segment to weld
            weld_ss = ss[ss[net.attr2].isin(all_ids)]
            weld = list(weld_ss.geometry)
            weld = _weld_MultiLineString(weld, skip_restr=net.skip_restr)
            
            # if the new segment if a LineString set the new, welded
            # geometry to the `keep_id` index of the dataframe
            if type(weld) == LineString:
                index = weld_ss.loc[(weld_ss[net.attr2] == keep_id)].index[0]
                gdf.loc[index, net.geo_col] = weld
                gdf.loc[index, phase] = phase_name
            
            # if the weld resulted in a MultiLineString remove ids from
            # from `drop_ids` and set to new for each n+1 new segment.
            if type(weld) == MultiLineString:
                unique_segs = len(weld)
                keeps_ids = [keep_id] + drop_ids[:unique_segs-1]
                index = list(weld_ss[weld_ss[net.attr2].isin(keeps_ids)].index)
                for idx, seg in enumerate(weld):
                    gdf.loc[index[idx], net.geo_col] = seg
                    gdf.loc[index, phase] = phase_name
                for idx in keeps_ids:
                    if idx in drop_ids:
                        drop_ids.remove(idx)
            
            # remove original segments used to create the new, welded
            # segment(s) from the full segments dataframe
            gdf = gdf[~gdf[net.attr2].isin(drop_ids)]
    
    gdf.reset_index(inplace=True, drop=True)
    if net.inter:
        gdf.to_file(net.inter+phase_name+net.file_type)
    
    return gdf


def check_cleanliness(net, df):
    """scan segments to determine cleanse_supercycle
    needs to be run again
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    df : geopandas.GeoDataFrame
        dataframe of geometries
    
    Returns
    -------
    cleaned : bool
        False (repeat another cycle) or True
        (all clean --> break out of loop)
    """
    
    # equal geometry check
    cleaned = streets_filter(net, df, equal_geom=True, scrub_check=True)
    if not cleaned:
        return cleaned
    
    # contained geometry check
    cleaned = streets_filter(net, df, contained_geom=True, scrub_check=True)
    
    return cleaned


def line_splitter(net, df, proj_init=None, save_keep=None,
                  inherit_attrs=False, calc_len=False, stage=None,
                  road_type='MTFCC'):
    """top level function for spliting line segements
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    df : geopandas.GeoDataFrame
        dataframe of line segments to split
    save_keep : str
        path to save the shapefile. Default is None.
    inherit_attrs : bool
     inherit attributes from the dominant line segment.
     Default is False.
    calc_len : bool
        calculated length and add column. Default is False.
    stage : str
        iteration of phase. Default is None.
    road_type : str
        column to use for grouping road types. Default is 'MTFCC'.
    
    Returns
    -------
    split_lines : geopandas.GeoDataFrame
        all line segments including unsplit lines
    """
    
    # it `net.mtfcc_split` is string put it into a list
    if not hasattr(net.mtfcc_split, '__iter__'):
        net.mtfcc_split = [net.mtfcc_split]
    
    # create subset of segments to split and not split
    if net.mtfcc_split_by and net.mtfcc_split:
        if type(net.mtfcc_split) == list:
            subset_codes = net.mtfcc_split_by + net.mtfcc_split
        elif type(net.mtfcc_split) == str:
            subset_codes = net.mtfcc_split_by + [net.mtfcc_split]
        non_subset = df[~df[road_type].isin(subset_codes)]
        df = df[df[road_type].isin(subset_codes)]
    
    if inherit_attrs:
        drop_cols = [net.len_col, net.geo_col, net.xyid, 'ring']
        attrs = [col for col in df.columns if not drop_cols.__contains__(col)]
        attr_vals = {attr:[] for attr in attrs}
    
    # Iterate over dataframe to find intersecting and split
    split_lines = []
    count_lines_split = 0
    for loi_idx in df.index:
        
        # Working with TIGER/Line *EDGES*
        if net.mtfcc_split_by and net.mtfcc_split:
            
            # if a line segment used for splitting
            # but not to be split itself
            if df[road_type][loi_idx] not in net.mtfcc_split:
                split_lines.extend([df[net.geo_col][loi_idx]])
                
                # fill dictionary with attribute values
                if inherit_attrs:
                    for attr in attrs:
                        fill_val = [df[attr][loi_idx]]
                        if attr == 'Phase':
                            fill_val = ['Split'+stage]
                        attr_vals[attr].extend(fill_val)
                continue
        
        # get segs from the dataset that intersect
        # with the Line Of Interest
        intersecting = get_intersecting_geoms(net, df1=df, geom1=loi_idx,
                                              wbool=False)
        intersecting = intersecting[intersecting.index != loi_idx]
        
        # Working with TIGER/Line *ROADS*
        if not net.mtfcc_split_by and not net.mtfcc_split:
            
            # if the LOI is an interstate only pass in
            # the ramps for splitting
            if df[road_type][loi_idx] == net.mtfcc_intrst:
                intersecting = intersecting[\
                                 (intersecting[road_type] == net.mtfcc_ramp)\
                               | (intersecting[road_type] == net.mtfcc_serv)\
                               | (intersecting[road_type] == net.mtfcc_intrst)]
            
            # if LOI not ramp and interstates in dataframe
            # filter them out
            elif list(intersecting[road_type]).__contains__(net.mtfcc_intrst)\
                and df[road_type][loi_idx] != net.mtfcc_ramp:
                intersecting = intersecting[intersecting[road_type]\
                                               != net.mtfcc_intrst]
        
        # if There are no intersecting segments
        if intersecting.shape[0] == 0:
            continue
        
        # ring road bool
        ring_road = literal_eval(df['ring'][loi_idx])
        
        # actual line split call happens here
        new_lines = _split_line(df[net.geo_col][loi_idx], loi_idx,
                                df=intersecting, ring_road=ring_road,
                                geo_col=net.geo_col)
        n_lines = len(new_lines)
        if n_lines > 1:
            count_lines_split += 1
        split_lines.extend(new_lines)
        
        # fill dictionary with attribute values
        if inherit_attrs:
            for attr in attrs:
                fill_val = [df[attr][loi_idx]]
                if attr == 'Phase' and n_lines != 1:
                    fill_val = ['Split'+stage]
                attr_vals[attr].extend(fill_val*n_lines)
    
    # create dataframe
    split_lines = gpd.GeoDataFrame(split_lines, columns=[net.geo_col])
    
    # fill dataframe with attribute values
    if inherit_attrs:
        for attr in attrs:
            split_lines[attr] = attr_vals[attr]
    if proj_init:
        split_lines = utils.set_crs(split_lines, proj_init=proj_init)
    if calc_len:
        split_lines = add_length(split_lines, len_col=net.len_col,
                                 geo_col=net.geo_col)
    
    # recombine EDGES subset and non subset segment lists
    if net.mtfcc_split_by and net.mtfcc_split:
        # combine newly split suset segments with all segments
        split_lines = split_lines.append(non_subset, sort=False)
        split_lines.reset_index(inplace=True, drop=True)
    split_lines = add_ids(split_lines, id_name=net.sid_name)
    
    # add updated xyid
    segm2xyid = utils.generate_xyid(df=split_lines, geom_type='segm',
                                    geo_col=net.geo_col)
    split_lines = utils.fill_frame(split_lines, col=net.xyid, data=segm2xyid)
    split_lines = label_rings(split_lines, geo_col=net.geo_col)
    
    # number of lines split
    net.lines_split = count_lines_split
    if save_keep:
        split_lines.to_file(save_keep+net.file_type)
    
    return split_lines


def _split_line(loi, idx, df=None, geo_col=None, ring_road=False):
    """middle level function for spliting line segements
    
    Parameters
    ----------
    loi : shapely.LineString
                line segment in question
    idx : int
                index number of the LOI
    df : geopandas.GeoDataFrame
                dataframe of line segments
    ring_road : bool
        (True) if ring road. (False) if not. Default is False.
    geo_col : str
        geometry column name. Default is None.
    
    Returns
    -------
    newLines : list
        list of new lines from one iteration of segment splitting
    """
    
    intersectinglines = df[df.index != idx] # all lines not LOI
    
    # Unary Union for intersection determination
    intersectinglines = intersectinglines[geo_col].unary_union
    
    # Intersections of LOI and the Unary Union
    breaks = loi.intersection(intersectinglines)
    
    # find and return points on the line to split if any exist
    unaltered, breaks,\
    ring_endpoint, basic_ring,\
    complex_ring = _find_break_locs(loi=loi, breaks=breaks,
                                    ring_road=ring_road)
    
    if unaltered:
        return unaltered
    
    # Line breaking
    if not type(breaks) == list:
        breaks = [breaks]
    new_lines = _create_split_lines(breaks=breaks, loi=loi,
                                    ring_road=ring_road, basic_ring=basic_ring,
                                    complex_ring=complex_ring,
                                    ring_endpoint=ring_endpoint)
    
    return new_lines


def _create_split_lines(breaks=None, ring_endpoint=None, loi=None,
                        ring_road=False, basic_ring=False, complex_ring=False):
    """deep function from splitting a single line segment
    along break points
    
    Parameters
    ----------
    breaks : list
        list of point to break a line along. Default is None.
    ring_endpoint : shapely.Point
        endpoint of a ring road. Default is None.
    loi : shapely.LineString
        line segment in question. Default is None.
    ring_road : bool
        is or is not a ring road. Default is False.
    basic_ring : bool
        is or is not a basic ring road. This indicates a 'normal' ring
        road in which there is one endpoint. Default is False.
    complex_ring : bool
        is or is not a complex ring road. This indicates a any situation
        not deemed a 'basic' ring. Default is False.
    
    Returns
    -------
    new_lines : list
        list of new lines generated from splitting
    """
    
    points = [Point(xy) for xy in breaks]
    
    # First coords of line
    coords = list(loi.coords)
    
    # Keep list coords where to cut (cuts = 1)
    cuts = [0] * len(coords)
    cuts[0] = 1
    cuts[-1] = 1
   
    # Add the coords from the points
    coords += [list(p.coords)[0] for p in points]
    cuts += [1] * len(points)
   
    # Calculate the distance along the line for each point
    dists = [loi.project(Point(p)) for p in coords]
    
    # sort the coords/cuts based on the distances
    # see http://stackoverflow.com/questions/6618515/
    #     sorting-list-based-on-values-from-another-list
    coords = [p for (d, p) in sorted(zip(dists, coords))]
    cuts = [p for (d, p) in sorted(zip(dists, cuts))]
    if ring_road: # ensure there is an endpoint for rings
        if basic_ring:
            coords = ring_endpoint + coords + ring_endpoint
        if complex_ring:
            archetelos = [loi.coords[0]] # beginning and ending of ring
            coords = archetelos + coords + archetelos
        cuts = [1] + cuts + [1]
    
    # generate the lines
    if cuts[-1] != 1: # ensure there is an endpoint for rings
        cuts += [1]
    new_lines = []
    for i in range(len(coords)-1):
        if cuts[i] == 1:
            # find next element in cuts == 1 starting from index i + 1
            j = cuts.index(1, i+1)
            new_line = LineString(coords[i:j+1])
            if new_line.is_valid:
                new_lines.append(new_line)
    return new_lines


def _find_break_locs(loi=None, breaks=None, ring_road=False,):
    """locate points along a line segment where breaks need to be made
    
    Parameters
    ----------
    loi : shapely.geometry.LineString
        Line of Interest
    breaks : list
        point to break a line. Default is None.
    ring_road : bool
        is or is not a ring road. Default is False.
    
    Returns
    -------
    unaltered : None or list
        list of one unaltered LineString
    breaks : None of list
        point to break a line. Default is None.
    ring_endpoint : shapely.Point
        endpoint of a ring road. Default is None.
    basic_ring : bool
        is or is not a basic ring road. This indicates a 'normal' ring
        road in which there is one endpoint.
    complex_ring : bool
        is or is not a complex ring road. This indicates a any
        situation not deemed a 'basic' ring.
    """
    
    intersection_type = type(breaks)
    unaltered = None
    ring_endpoint = None
    basic_ring = False
    complex_ring = False
    
    # Case 1
    # - Single point from a line intersects the LOI
    # loop roads & 'typical' intersections
    if intersection_type == Point:
        ring_endpoint = [breaks]
        if ring_road == False:
            if breaks == loi.boundary[0]\
            or breaks == loi.boundary[1]:
                return [loi], None, None, None, None
        if ring_road:
            basic_ring = True
           # Do nothing, return the exact ring geometry
            return [loi], None, None, None, None
        else:
            breaks = _make_break_locs(breaks=breaks, standard=True)
    
    # Case 2 
    # - Multiple points from one line intersect the LOI
    # horseshoe roads, multiple intersections of one line and LOI
    elif intersection_type == MultiPoint:
        if ring_road:
            complex_ring = True
            breaks = _make_break_locs(breaks=breaks)
        # horseshoe
        elif breaks == loi.boundary:
            return [loi], None, None, None, None
        else:
            breaks = _make_break_locs(loi=loi, breaks=breaks)
    
    # Case 3
    # - Overlapping line segments along one stretch of road
    # multiple names, etc. for a section of roadway which was then
    # digitized as separate, stacked entities.
    elif intersection_type == LineString:
        breaks = _make_break_locs(loi=loi, breaks=breaks, line=True)
    
    # Case 4
    # - Overlapping line segments along multiple stretches of
    # road multiple names, etc. for multiple sections of roadway which
    # were then digitized as separate, stacked entities.
    elif intersection_type == MultiLineString:
        breaks = _make_break_locs(loi=loi, breaks=breaks, mline=True)
    
    # Case 5
    # - Complex intersection of points and Lines
    # anomaly in digitization / unclear
    elif intersection_type == GeometryCollection:
        # points and line in the geometry collection
        pts_in_gc = []
        lns_in_gc = []
        
        #if only one line intersection with a point intersection
        multiple_line_intersections = False
        for idx, geom in enumerate(breaks):
            # collect points
            if type(geom) == Point or type(geom) == MultiPoint:
                pts_in_gc.append(geom)
            # collect line(s)
            if type(geom) == LineString or type(geom) == MultiLineString:
                lns_in_gc.append(geom)
        
        #get split indices in line based on touching geometry
        split_index = [] 
        iter_limit = len(lns_in_gc)-1
        for i in range(iter_limit):
            j = i+1
            current_geom = lns_in_gc[i]
            next_geom = lns_in_gc[j]
            # comparing incremental geometry pairs
            # if touching: do nothing
            if current_geom.touches(next_geom):
                continue
            else: # if don't touch add a split index
                split_index.append(j)
                multiple_line_intersections = True
        
        # if there are multiple line intersections between two
        # lines split the segments at the line intersections
        if multiple_line_intersections:
            split_list = []
            for split in split_index:
                if split == split_index[0]:
                    prev_split = split
                    # first split
                    section = lns_in_gc[:split]
                else:
                     # 2nd to n-1 split
                    section = lns_in_gc[prev_split:split]
                # add split line segment
                split_list.append(section)
                # only one split
                if split_index[0] == split_index[-1]:
                    split_list.append(lns_in_gc[split:])
                 # last split
                elif split == split_index[-1]:
                    split_list.append(lns_in_gc[split:])
            lns_in_gc = split_list
        
        # otherwise if there are not multiple line intersections...
        if not multiple_line_intersections:
            welded_line = _weld_MultiLineString(lns_in_gc)
            pts_in_gc.extend([welded_line.boundary[0],
                              welded_line.boundary[1]])
        elif multiple_line_intersections:
            for geoms in lns_in_gc:
                if len(geoms) == 1:
                    pts_in_gc.extend([geoms[0].boundary[0],\
                                      geoms[0].boundary[1]])
                else:
                    welded_line = _weld_MultiLineString(geoms)
                    pts_in_gc.extend([welded_line.boundary[0],
                                      welded_line.boundary[1]])
        
        breaks = _make_break_locs(loi=loi, breaks=pts_in_gc)
    
    return unaltered, breaks, ring_endpoint, basic_ring, complex_ring


def _make_break_locs(breaks=None, standard=False, loi=None,
                     line=False, mline=False):
    """record the points along a line where breaks needs to be
    made after finding them
    
    Parameters
    ----------
    breaks : shapely.Point or shapely.LineString
        object to use for breaking the segment. Default is None.
    standard : bool
        this indicates a single point break. Default is False.
    loi : shapely.LineString coordinates
        coordinates along a line. Default is None.
    line : bool
        is a LineString. Default is False.
    mline : bool
        is a MultiLineString. Default is False.
    
    Returns
    -------
    break_points : list
        geometries of points to break a line
    """
    
    if breaks and standard:
        break_points = [Point(breaks.coords[:][0])]
    
    elif breaks and not standard:
        if line:
            breaks = [breaks.boundary[0], breaks.boundary[1]]
        if mline:
            lns_in_mls = [l for l in breaks]
            
            # created welded line, but only referencing
            # for break location
            welded_line = _weld_MultiLineString(lns_in_mls)
            breaks = [welded_line.boundary[0], welded_line.boundary[1]]
        break_points = [Point(point.coords[:][0]) for point in breaks]
    
    return break_points


def _weld_MultiLineString(multilinestring, weld_multi=True, skip_restr=True):
    """weld a shapely.MultiLineString into a shapely.LineString
    
    Parameters
    ----------
    multilinestring : shapely.geometry.MultiLineString
        segment (collection) to weld.
    weld_multi :bool
        if welded line is still a multiline segment, then determine if
        the segments of the multiline are almost equal. Default is True.
    skip_restr : bool
        skip re-welding restricted segments. Default is False.
    
    Returns
    -------
    welded : shapely.geometry.LineString
        freshly welded segment (collection).
    """
    
    welded = linemerge(multilinestring)
    
    # Due to minute rounding (.00000001 meters) some line vertices can
    # be off thus creating a MultiLineString where shapely thinks two
    # LineString objects don't actually touch where, in fact, they do.
    # The following loop iterates through each pair of LineString
    # objects sequentially to  determine if their endpoints are
    # 'almost equal' instead of exactly equal. When the endpoints are
    # 'almost equal' the starting point of the second line is duplicated
    # as the ending point of the first line before the lines are
    # welded together.
    if type(welded) == MultiLineString and weld_multi and skip_restr:
        line_count = len(welded)
        new_lines = {}
        for line1 in range(line_count):
            for line2 in range(line1+1, line_count):
                
                # geometries
                L1, L2 = welded[line1], welded[line2]
                # starting and endpoints
                sp1, ep1 = L1.boundary[0], L1.boundary[1]
                sp2, ep2 = L2.boundary[0], L2.boundary[1]
                
                # if equal move along
                if ep1.equals(sp2) or sp1.equals(ep2):
                    continue
                
                # if either sets are almost equal pass along the
                # altered first line and the original second line
                if ep1.almost_equals(sp2) or sp1.almost_equals(ep2):
                    if ep1.almost_equals(sp2):
                        new_line = LineString(L1.coords[:-1] + L2.coords[:1])
                    if sp1.almost_equals(ep2):
                        new_line = LineString(L2.coords[-1:] + L1.coords[1:])
                    new_lines[line1] = new_line
        
        # convert welded multiline to list
        welded = list(welded)
        for idx, line in list(new_lines.items()):
            welded[idx] = line
        
        # re-weld
        welded = linemerge(welded)
    
    return welded


def _drop_geoms(net, gdf, geoms, series=False,
                save_kept=None, save_dropped=None):
    """Drop a subset of geometries from a geopandas dataframe
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    gdf : geopandas.GeoDataFrame
        dataframe of geometries to search
    geoms : list
        either a list or a list of lists
    series : bool
        search a geoseries. Default is False.
    save_kept : str
        path to save filtered geometries. Default is None.
    save_dropped : str
        path to save dropped geometries. Default is None.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        kept geometries in a dataframe
    """
    
    if series:
        drop_geoms = geoms
    else:
        drop_geoms = set([item for sublist in geoms for item in sublist])
    
    if None in drop_geoms:
        drop_geoms.remove(None)
    dropped_geoms = gdf[gdf.index.isin(drop_geoms)]
    
    if save_dropped:
        dropped_geoms.to_file(save_dropped+file_type)
    gdf = gdf[~gdf.index.isin(drop_geoms)]
    
    if save_kept:
        gdf.to_file(save_kept+net.file_type)
    
    return gdf


def create_node(x,y):
    """create a node along the network
    
    Parameters
    ----------
    x : float or int
        x coordinate of point
    y : loat or int
        y coordinate of point
    
    Returns
    -------
    instantiated shapely.Point
    """
    
    return Point(list(zip(x, y))[0])


def create_edges(s2n, ndf, c1, c2):
    """create a edge along the network
    
    Parameters
    ----------
    s2n : list
        segment2node relationship list
    ndf : geopandas.GeoDataFrame
        nodes dataframe
    c1 : str
        column 1 name
    c2 : str
        column 2 name
    
    Returns
    -------
    s2e : dict
        segment2edge dictionary
    """
    
    s2e = {}
    
    for seg, nodes in list(s2n.items()):
        incident_vertices = []
        for node in nodes:
            incident_vertices.append(ndf.loc[(ndf[c1] == node), c2].values[0])
        s2e[seg] = LineString(incident_vertices)
    
    return s2e


def seg_welder(net, frame):
    """return a copy of the welded line segment within a frame
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    frame : geopandas.GeoDataframe
        frame of line segments
    
    Returns
    -------
    f : geopandas.GeoDataframe
        frame of line segments
    """
    
    f = frame.copy()
    welded_segs = 0
    
    for idx in f.index:
        if type(f[net.geo_col][idx]) == MultiLineString:
            welded_segs += 1
            f[net.geo_col][idx] = _weld_MultiLineString(f[net.geo_col][idx])
    
    del frame
    
    net.welded_mls = welded_segs
    
    return f


def _seg_partition(lst, n):
    """Partition one road segment into (n+1) smaller segments with (n)
    
    Parameters
    ----------
    lst : list
        list coordinates in LineString
    n : int
        number of break points
    
    Returns
    -------
    node_overlap : list
        overlapping endpoints of segments
    """
    
    if len(lst) < n:
        n = len(lst) 
    div = len(lst) / float(n)
    xoverlap = [lst[int(round(div*i)):int(round(div*(i+1)))] for i in range(n)]
    node_overlap = []
    
    for idx, part in enumerate(xoverlap):
        if idx == len(xoverlap)-1:
            node_overlap.append(part)
        else:
            p1, p2 = part, [xoverlap[idx+1][0]]
            node_overlap.append(p1+p2)
    
    node_overlap = [section[-1] for section in node_overlap]
    
    return node_overlap


def series_search(net, df, seg=None, scanned=None, equal_geom=False,
                  contained_geom=False, covered_geom=False, drop_over=False,
                  scrub_check=False):
    """compare a single line segment to a GeoSeries
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    df : geopandas.GeoDataFrame
        dataframe of geometries to search
    seg : int
         geometry id in the dataframe. Default is None.
    equal_geom : bool
        compare for exactly equal geometries. Default is False.
    contained_geom : bool
        compare for contained geometries. Default is False.
    covered_geom : bool
        compare for partially covered geometries. Default is False.
    drop_over : bool
        drop covering geometries. Default is False.
    scanned : set
        already scanned geometries. Default is None.
    scrub_check : bool
        performing a scrub check (True), or not (False).
        Default is False.
    
    Returns
    -------
    drop_seg : int
        segment to drop/mask out
    """
    
    drop_segs = [None]
    
    # intersecting bool and intersecting segments subset
    i_bool, i_segs = get_intersecting_geoms(net, df1=df, geom1=seg)
    
    if i_segs.shape[0] > 1:
        if equal_geom:
            i_segs = i_segs[i_segs.index != seg]
            drop_segs = [None]
            
            for i_idx in i_segs.index:
                if i_segs[net.geo_col][i_idx].equals(df[net.geo_col][seg])\
                and i_idx not in scanned:
                    
                    # use for determining whether to run
                    # another iteration of supercleanse
                    if scrub_check:
                        return False
                    
                    # same geometry; different MTFCC value;
                    # s1 is interstate;
                    if df[net.attr1][seg] != i_segs[net.attr1][i_idx]\
                    and df[net.attr1][seg] == net.mtfcc_split\
                    and i_segs[net.attr1][i_idx] != net.mtfcc_split:
                        drop_segs.append(i_idx)
                    
                    # same geometry; different MTFCC value;
                    # s2 is interstate;
                    if df[net.attr1][seg] != i_segs[net.attr1][i_idx]\
                    and df[net.attr1][seg] != net.mtfcc_split\
                    and i_segs[net.attr1][i_idx] == net.mtfcc_split:
                        drop_segs.append(seg)
                    
                    # same geometry; different MTFCC value;
                    # neither s1 or s2 is interstate
                    if df[net.attr1][seg] != i_segs[net.attr1][i_idx]\
                    and df[net.attr1][seg] != net.mtfcc_split\
                    and i_segs[net.attr1][i_idx] != net.mtfcc_split:
                        drop_segs.append(i_idx)
                    
                    # same geometry same MTFCC value
                    if df[net.attr1][seg] == i_segs[net.attr1][i_idx]\
                    and df[net.attr1rank][seg]\
                    == i_segs[net.attr1rank][i_idx]:
                        drop_segs.append(i_idx)
                    
                    # same geometry; different MTFCC value;
                    # not interstates
                    if df[net.attr1][seg] != i_segs[attr1][i_idx]\
                    and df[net.attr1][seg] != net.mtfcc_split\
                    and i_segs[net.attr1][i_idx] != net.mtfcc_split:
                        
                        if net.attr1 == 'MTFCC':
                            # keep the less likely to have
                            # a higher speed limit
                            mtfcc1, mtfcc2 = int(df[net.attr1][seg][1:]),\
                                             int(i_segs[net.attr1][i_idx][1:])
                            if df[net.attr1rank][seg]\
                            > i_segs[net.attr1rank][i_idx]:
                                drop_segs.append(seg)
                            else:
                                drop_segs.append(i_idx)
        
        if contained_geom:
            i_segs = i_segs[i_segs.index != seg]
            drop_segs = [None]
            for i_idx in i_segs.index:
                
                # seg1 contains seg2 and both have same attribute value
                if df[net.geo_col][seg].contains(i_segs[net.geo_col][i_idx])\
                and df[net.attr1rank][seg] == i_segs[net.attr1rank][i_idx]\
                and i_idx not in scanned:
                    if scrub_check:
                        return False
                    if drop_over:
                        drop_segs.append(seg)
                    else:
                        drop_segs.append(i_idx)
                
                # seg2 contains seg1 and both have same attribute value
                if i_segs[net.geo_col][i_idx].contains(df[net.geo_col][seg])\
                and i_segs[net.attr1rank][i_idx] == df[net.attr1rank][seg]\
                and not scanned.__contains__(i_idx):
                    if scrub_check:
                        return False
                    if drop_over:
                        drop_segs.append(i_idx)
                    else:
                        drop_segs.append(seg)
                
                # seg1 contains seg2 but have different attribute value
                if df[net.geo_col][seg].contains(i_segs[net.geo_col][i_idx])\
                and df[net.attr1rank][seg] != i_segs[net.attr1rank][i_idx]\
                and not scanned.__contains__(i_idx):
                    if scrub_check:
                        return False
                    if drop_over:
                        drop_segs.append(seg)
                    else:
                        drop_segs.append(i_idx)
                
                # seg2 contains seg1 but have different attribute value
                if i_segs[net.geo_col][i_idx].contains(df[net.geo_col][seg])\
                and i_segs[net.attr1rank][i_idx] == df[net.attr1rank][seg]\
                and not scanned.__contains__(i_idx):
                    if scrub_check:
                        return False
                    if drop_over:
                        drop_segs.append(i_idx)
                    else:
                        drop_segs.append(seg)
        
        if covered_geom:
            i_segs = i_segs[i_segs.index != seg]
            drop_segs = [None]
            
            for i_idx in i_segs.index:
                # same geometry and same attribute value
                isec = df[net.geo_col][seg]\
                       .intersection(i_segs[net.geo_col][i_idx])
                inter_type = type(isec)
                equal_attr = df[net.attr1rank][seg]\
                             == i_segs[net.attr1rank][i_idx]
                diff_attr = df[net.attr1rank][seg]\
                            != i_segs[net.attr1rank][i_idx]
                
                if (inter_type == MultiLineString and equal_attr)\
                or (inter_type == GeometryCollection\
                and len(isec) > 0 and equal_attr):
                    if scrub_check:
                        return False
                    drop_segs.extend([seg, i_idx])
                
                # same geometry and different attribute value
                if (inter_type == MultiLineString and diff_attr)\
                or (inter_type == GeometryCollection\
                and len(isec) > 0 and diff_attr):
                    if scrub_check:
                        return False
                    drop_segs.extend([seg, i_idx])
    
    if scrub_check:
        return True
    
    return drop_segs


def drop_micro_slivers(gdf, tol=0.01, len_col=None, id_name=None):
    """Drop micro sliver that are gnerated following line splits
    
    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        dataframe of line segments
    tol : float or int
        drop all segment with a length less than this. Default is 0.01.
    len_col : str
        length column name in dataframe. Default is None.
    id_name  : str
        id column name in dataframe. Default is None.
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        updated dataframe of line segments 
    """
    
    predrop_shape = gdf.shape
    gdf = gdf[gdf[len_col] > tol]
    postdrop_shape = gdf.shape
    slivers = predrop_shape[0] - postdrop_shape[0]
    gdf.reset_index(drop=True, inplace=True)
    gdf = add_ids(gdf, id_name=id_name)
    
    return gdf


def streets_filter(net, df, series=False, equal_geom=False,
                   contained_geom=False, covered_geom=False,
                   drop_over=False, save_keep=None,
                   save_drop=None, scrub_check=False):
    """top level function to filter and pare down street segments.

    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    df : geopandas.GeoDataFrame
        streets dataframe
    series  : bool
        search a geoseries. Default is False.
    equal_geom : bool
        drop all but one of equal geometries. Default is False.
    contained_geom : bool
        drop all contained geometries. Default is False.
    covered_geom : bool
        drop all covered geometries. Default is False.
    drop_over : bool
        drop covering geometries. Default is False.
    save_keep : str
        path to save filtered geometries. Default is None.
    save_drop : str
        path to save dropped geometries. Default is None.
    scrub_check : bool
        performing a scrub check (True), or not (False).
        Default is False.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        dataframe of filtered streets
    """
    
    unclean = True
    scanned = set()
    drop = set()
    
    for segm in df.index:
        if segm not in scanned:
            d = series_search(net, df, seg=segm, scanned=scanned,
                              equal_geom=equal_geom, drop_over=drop_over,
                              contained_geom=contained_geom,
                              covered_geom=covered_geom,
                              scrub_check=scrub_check)
            
            if scrub_check:
                cleaned = d
                if not cleaned:
                    return cleaned
                else:
                    continue
            
            if type(d) == list:
                drop.update(d)
            
            else: 
                drop.add(d)
            scanned.add(segm)
    
    if scrub_check:
        return cleaned
    
    df = _drop_geoms(net, df, drop, series=series,
                     save_kept=save_keep, save_dropped=save_drop)
    
    return df



###############################################################################
################ spaghetti.SpaghettiNetwork functionality #####################
###############################################################################



def extract_nodes(net):
    """extract n_ids from line segments and return
    them in a geodataframe
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    
    Returns
    -------
    nodedf : geopandas.GeoDataFrame
        node dataframe
    """
    
    def _drop_covered_nodes(net, nodedf):
        """keep only the top node in stack of overlapping nodes
        
        Parameters
        ----------
        net : spaghetti.SpaghettiNetwork
        nodedf : geopandas.GeoDataFrame
            node dataframe
        
        Returns
        -------
        nodedf : geopandas.GeoDataFrame
            updated node dataframe
        """
        
        scanned = set()
        drop = set()
        
        # keep only the top node in a node stack.
        for n in nodedf.index:
            if n not in scanned:
                unscanned = nodedf[~nodedf.index.isin(scanned)]
                xyid = unscanned[net.xyid][n]
                
                # intersecting node idxs
                iidxs = list(unscanned[unscanned[net.xyid] == xyid].index)
                iidxs.remove(n)
                if iidxs:
                    drop.update(iidxs), scanned.update(drop)
                scanned.add(n)
        nodedf = _drop_geoms(net, nodedf, drop, series=True)
        nodedf = add_ids(nodedf, id_name=net.nid_name)
        
        return nodedf
    
    sdf = net.s_data
    nodes = []
    
    # create n_ids and give the segment attribute data
    for seg in sdf.index:
        if sdf['ring'][seg] == 'True':
            xs, ys = sdf[net.geo_col][seg].coords.xy
            nodes.append(create_node(xs, ys))
        else:
            nodes.extend([sdf[net.geo_col][seg].boundary[0],
                          sdf[net.geo_col][seg].boundary[1]])
    nodedf = gpd.GeoDataFrame(nodes, columns=[net.geo_col])
    nodedf = add_ids(nodedf, id_name=net.nid_name)
    
    if net.proj_trans:
        nodedf = utils.set_crs(nodedf, proj_init=net.proj_trans)
    
    # Give an initial string 'xy' ID 
    prelim_xy_id = utils.generate_xyid(df=nodedf, geom_type='node',
                                       geo_col=net.geo_col)
    nodedf = utils.fill_frame(nodedf, idx=net.nid_name,
                              col=net.xyid, data=prelim_xy_id)
    
    # drop all node but the top in the stack
    nodedf = _drop_covered_nodes(net, nodedf)
    nodedf.reset_index(drop=True, inplace=True)
    
    return nodedf


def associate(primary=None, secondary=None, assoc=None,
              initial_weld=False, net=None, df=None, ss=None):
    """create 2 dicts of neighor relationships (x2y and y2x). *OR*
    create one list of x2y neighor relationships.
    
    Parameters
    ----------
    primary : list
        primary data in the form - [idx, [xyID1, xyID2,...]].
        Default is None.
    secondary : list
        secondary data in the form [idx, [xyID1, xyID2,...]].
        Default is None.
    assoc : str
        either node2seg or seg2node. Default is None.
    initial_weld : bool
        welding subset of restricted access road segments. Used in
        cleanse_supercycle(). Default is False.
    net : spaghetti.SpaghettiNetwork
        Default is None.
    df : geopandas.GeoDataFrame
        restricted streets susbet dataframe. Default is None.
    ss : geopandas.GeoDataFrame
        subset of restricted streets susbet. Default is None.
    
    Returns
    -------
    segm_dict : dict
        neighoring elements in the form - {seg1, [node1, node2]}
    node_dict : dict
        neighoring elements in the form - {node1, [seg1, seg2]}
    topos_list : list
        neighoring elements in the form - [x1, [y1, y2]]
    """
    
    if initial_weld:
        segm_dict = {}
        
        for idx in df.index:
            neigh = [ss[net.tnidf][idx], ss[net.tnidt][idx]]
            segm_dict[ss[net.attr2][idx]] = neigh
        
        # get nodes
        ss_nodes = set()
        
        for sidx, nidx in list(segm_dict.items()):
            for nx in nidx:
                ss_nodes.add(nx)
        node_dict = {}
        
        for node_idx in ss_nodes:
            node_dict[node_idx] = set()
            for seg_idx, nodes_idx in list(segm_dict.items()):
                if node_idx in nodes_idx:
                    node_dict[node_idx].add(seg_idx)
        
        return segm_dict, node_dict
    
    topos_list = []
    
    for primary_idx, primary_info in enumerate(primary):
        topos = [primary_idx, []]
        
        for secondary_idx, secondary_info in enumerate(secondary):
            secondary_idxs = []
            
            # first and last point of the segment in string format
            # for primary_info in 'segm2node' and
            # secondary_info in 'node2segm'
            if assoc == 'segm2node':
                if secondary_info[1][0] == primary_info[1][0]\
                or secondary_info[1][0] == primary_info[1][-1]:
                    secondary_idxs.append(secondary_idx)
            
            if assoc == 'node2segm':
                if primary_info[1][0] == secondary_info[1][0]\
                or primary_info[1][0] == secondary_info[1][-1]:
                    secondary_idxs.append(secondary_idx)
            
            topos[1].extend(secondary_idxs)
        
        topos_list.extend([topos])
    
    return topos_list


def get_neighbors(x2y, y2x, astype=None):
    """get all neighboring graph elements of the same type
    
    Parameters
    ----------
    x2y : list or dict
        element type1 to element type2 crosswalk
    y2x : list or dict
        element type2 to element type1 crosswalk
    astype : list or dict
        return the lookup as either type. Default is None.
    Returns
    -------
    x2x : list or dict
        element type1 to element type1 crosswalk OR element type2 to
        element  type2 crosswalk in the form
        - {x1, [x2,x3]} *OR* [x1, [x2,x3]]
    """
    
    if not astype:
        raise Exception('`astype` parameter must be set.')
    
    elif astype == dict:
        x2x = {}
        for k, vn in list(x2y.items()):
            x2x[k] = set()
            for v in vn:
                x2x[k].update(y2x[v])
                x2x[k].discard(k)
    
    elif astype == list:
        x2x = []
        for (k, vn) in x2y:
            x2x.append([k, set()])
            for v in vn:
                x2x[k][1].update(y2x[v][1])
                x2x[k][1].discard(k)
        x2x = [[k, list(v)] for (k,v) in x2x]
    
    else:
        raise Exception(str(type), 'not a valid type for `astype` parameter.')
    
    return x2x


def xwalk(df, c1=None, c2=None, stipulation=None, geo_col=None):
    """create adjacency cross walks as lists
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        geometry dataframe
    c1 : str
        column name. Default is None.
    c2 : str
        column name. Default is None.
    stipulation : str
        Default is None.
    geo_col : str
        name of geometry column.
    
    Returns
    -------
    xw : list
    """
    
    if c2 == 'nodeNeighs' or c2 == 'segmNeighs':
        xw = [[df[c1][ix],literal_eval(df[c2][ix])] for ix in df.index]
    
    if c2 == 'degree' or c2 == 'length' or c2 == 'TLID':
        xw = [[df[c1][ix],df[c2][ix]] for ix in df.index]
    
    if c2 == geo_col and not stipulation:
        xw = [[df[c1][ix],df[c2][ix]] for ix in df.index]
    
    if c2 == geo_col and stipulation == 'coords':
        xw = [[df[c1][ix],df[c2][ix].coords[:]] for ix in df.index]
    
    return xw


def assert_2_neighs(net):
    """
    1. raise an error if a segment has more that 2 neighbor nodes.
    2. if the road has one neighbor node then it is a ring road. In 
    this case give the ring road a copy of the one nighbor.
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    """
    
    more_2_neighs = [k for (k,v) in net.segm2node if len(v) > 2]
    
    if more_2_neighs:
        raise ValueError('Adjacency value corruption. The segments '\
                         + 'listed below are incident with more than '\
                         + 'two nodes.\n\n'\
                         + 'Problem segment IDs: ' + str(more_2_neighs))
    
    rings = [k for (k,v) in net.segm2node if len(v) < 2]
    
    if rings:
        for ring in rings:
            n1 = net.segm2node[ring][1][0]
            net.segm2node[ring][1] = [n1,n1]
    
    return net


def get_cc_len(net, len_col=None):
    """return the geodataframe with the length of each associated
    connected component in a new column.
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    len_col : str
        name of length column. Default is None.
    
    Returns
    -------
    cc_lens : dict
        {ID:length} for each connected component in graph
    """
    
    net.s_data['ccLength'] = np.nan
    cc_lens = {}
    
    for (k,v) in net.segm_cc:
        new_v, segment_ids = v, v
        new_v = net.s_data[net.s_data[net.sid_name].isin(new_v)]
        new_v = new_v[len_col].sum()
        net.s_data.loc[net.s_data[net.sid_name].isin(v),'ccLength'] = new_v
        cc_lens[k] = [new_v, segment_ids]
    
    return cc_lens


def get_largest_cc(ccs, smallkeys=True):
    """connected components object
    
    Parameters
    ----------
    ccs : list
        list of connected components
    smallkeys : bool
        Default is True.
    
    Returns
    -------
    [largestKey, largestValues]: largest lookup
    nonLargest : list
        all non largest keys
    """
    
    largest = max(ccs, key=lambda k: len(k[1]))
    largest_key = largest[0]
    largest_values = largest[1]
    
    if smallkeys:
        
        non_largest = []
        for (cck, ccvs) in ccs:
            if cck is not largest_key:
                non_largest.extend(ccvs)
        
        return [largest_key, largest_values], non_largest
    
    return [largest_key, largest_values]


def update_adj(net, seg_keys, node_keys):
    """update adjacency relationships between segments and nodes.
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    seg_keys : list
        segment keys to remove from adjacency
    node_keys : list
        node keys to remove from adjacency
    """
    
    # update all crosswalk dictionaries
    net.segm2segm = remove_adj(net.segm2segm, seg_keys)
    net.segm2node = remove_adj(net.segm2node, seg_keys)
    net.node2node = remove_adj(net.node2node, node_keys)
    net.node2segm = remove_adj(net.node2segm, node_keys)
    
    # Keep only the largest connected components
    net.segm_cc = net.largest_segm_cc
    net.node_cc = net.largest_node_cc
    
    # Set component ID to dataframe
    net.s_data = net.s_data[net.s_data[net.sid_name].isin(net.segm_cc[1])]
    net.s_data.reset_index(drop=True, inplace=True)
    net.n_data = net.n_data[net.n_data[net.nid_name].isin(net.node_cc[1])]
    net.n_data.reset_index(drop=True, inplace=True)


def remove_adj(x2x, remove):
    """remove adjacent elements from list of ids
    
    Parameters
    ----------
    x2x : list
        x2x relationship list
    remove : list
        keys to remove from list
    
    Returns
    -------
    x2x : list
        updated x2s relationship list
    """
    
    return [[k, vs] for (k,vs) in x2x if k not in set(remove)]


def geom_assoc(net):
    """Associate nodes and segments with geometry
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    """
    
    net.segm2geom = xwalk(net.s_data, c1=net.sid_name,
                          c2=net.geo_col, geo_col=net.geo_col)
    net.node2geom = xwalk(net.n_data, c1=net.nid_name,
                          c2=net.geo_col, geo_col=net.geo_col)


def coords_assoc(net):
    """Associate nodes and segments with coordinates
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    """
    
    net.segm2coords = xwalk(net.s_data, c1=net.sid_name, c2=net.geo_col,
                            geo_col=net.geo_col, stipulation='coords')
    net.node2coords = xwalk(net.n_data, c1=net.nid_name, c2=net.geo_col,
                            geo_col=net.geo_col, stipulation='coords')


def euc_calc(net, col=None):
    """Calculate the euclidean distance between two line endpoints
    for each line in a set of line segments
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    col : str
        new column name. Default is None.
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        updated segments dataframe
    """
    
    net.s_data[col] = np.nan
    for (seg_k, (n1,n2)) in net.segm2node:
        p1, p2 = net.node2coords[n1][1][0], net.node2coords[n2][1][0]
        ed = _euc_dist(p1, p2)
        net.s_data.loc[(net.s_data[net.sid_name] == seg_k), col] = ed
    
    return net.s_data


def calc_valency(net, col=None):
    """calculate the valency fo each node and return a lookup
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    col : str
        node neighbors. Default is None.
    
    Returns
    -------
    n2d : list
        node2degree lookup
    """
    
    n2d = []
    for (node, segs) in net.node2segm:
        loops = 0
        for s in segs:
            neighs = literal_eval(\
                     net.s_data.loc[(\
                     net.s_data[\
                     net.sid_name] == s), col].values[0])
            if neighs[0] != neighs[1]:
                continue
            if neighs[0] == neighs[1]:
                loops += 1
        degree = len(segs) + loops
        n2d.append([node, [degree]])
    
    return n2d


def branch_or_leaf(net, geom_type=None):
    """Define each graph element (either segment or node) as either
    branch or leaf. Branches are nodes with degree 2 or higher, or
    segments with both incident nodes of degree 2 or higher
    (a.k.a. internal elements). Leaves are nodes with degree 1 or
    less, or segments with one incident node of degree 1 (a.k.a 
    external elements). Branches are 'core' elements, while leaves
    can be thought of as 'dead-ends'.
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    geom_type : str
        segm or node
    
    Returns
    -------
    geom2ge : list
        geometry idx to graph element type crosswalk
    """
    
    if geom_type == 'segm':
        id_list = net.segms_id_list
    else:
        id_list = net.n_ids
    
    geom2ge = []
    for idx in id_list:
        if geom_type == 'segm':
            n1 = net.segm2node[idx][1][0]
            n2 = net.segm2node[idx][1][1]
            n1d = net.node2degree[n1][1][0]
            n2d = net.node2degree[n2][1][0]
            
            if n1d == 1 or n2d == 1:
                graph_element = 'leaf'
            else:
                graph_element = 'branch'
        if geom_type == 'node':
            nd = net.node2degree[idx][1][0]
            if nd == 1:
                graph_element = 'leaf'
            else:
                graph_element = 'branch'
        geom2ge.append([idx, graph_element])
    
    return geom2ge


def simplify(net):
    """remove all non articulation objects
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    
    Returns
    -------
    segs : geopandas.GeoDataFrame
        simplified segments dataframe
    """
    
    # locate all non-articulation points
    na_objs = _locate_naps(net)
    
    # remove all non-articulation points
    segs = _simplifysegs(net, na_objs)
    
    return segs


def get_stats_frame(net):
    """create dataframe with descriptive network
    statistics as an attribute
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    """
    
    
    def _get_len_info(net, stat):
        """
        Parameters
        ----------
        stat : str
            either 'min' or 'max'
        
        Returns
        -------
        info : str
            segment id and length
        """
        if stat == 'min':
            lstat = net.s_data[net.len_col].min()
        
        else:
            lstat = net.s_data[net.len_col].max()
        stat_loc = (net.s_data[net.len_col] == lstat)
        sid_vals = net.s_data.loc[stat_loc, net.sid_name].values[:]
        info = str(sid_vals) + ' - ' + str(lstat)
        
        return info
    
    
    from collections import OrderedDict
    network_stats = ['Type', 'Measure', 'Value']
    
    # Node stats
    stats_dict = OrderedDict()
    stats_dict['n Nodes'] = net.n_node, 'Node'
    stats_dict['Min Node Degree'] = net.min_node_degree, 'Node'
    stats_dict['Max Node Degree'] = net.max_node_degree, 'Node'
    stats_dict['Mean Node Degree'] = net.mean_node_degree, 'Node'
    stats_dict['StD Node Degree'] = net.std_node_degree, 'Node'
    
    # Segment stats
    stats_dict['n Segments'] = net.n_segm, 'Segment'
    stats_dict['Total Length'] = net.network_length, 'Segment'
    stats_dict['Min Length'] = _get_len_info(net, 'min'), 'Segment'
    stats_dict['Max Length'] = _get_len_info(net, 'max'), 'Segment'
    stats_dict['Mean Length'] =net.s_data[net.len_col].mean(), 'Segment'
    stats_dict['StD Length'] =net.s_data[net.len_col].std(), 'Segment'
    stats_dict['Radius'] = str(net.radius[0])+' - '\
                           +str(net.radius[1]), 'Segment'
    stats_dict['Diameter'] = str(net.diameter[0])+' - '\
                             +str(net.diameter[1]), 'Segment'
    stats_dict['Min Sinuosity'] = net.min_sinuosity, 'Segment'
    stats_dict['Max Sinuosity'] = net.max_sinuosity, 'Segment'
    stats_dict['Mean Sinuosity'] = net.network_mean_sinuosity, 'Segment'
    stats_dict['StD Sinuosity'] = net.network_std_sinuosity, 'Segment'
    stats_dict['Alpha'] = net.alpha, 'Segment'
    stats_dict['Beta'] = net.beta, 'Segment'
    stats_dict['Eta'] = net.eta, 'Segment'
    stats_dict['Gamma'] = net.gamma, 'Segment'
    stats_dict['Circuity'] = net.circuity, 'Segment'
    stats_dict['MTFCC Entropy'] = net.entropy_mtfcc, 'Segment'
    
    # Object size
    stats_dict['Object Size (GB)'] = net.actual_total_size, 'Object'
    
    # set empty shell to fill
    empty_shell = np.empty((len(stats_dict), len(network_stats)))
    
    # dataframe containing descriptive network statistics
    net.network_stats = gpd.GeoDataFrame(empty_shell, columns=network_stats)
    
    # fill dataframe record by record
    for idx, (stat_name, (stat_value, stat_type))\
    in enumerate(stats_dict.items()):
        net.network_stats.iloc[idx] = stat_type, stat_name, stat_value
    
    # set dataframe index as type
    net.network_stats = net.network_stats.set_index(['Type'])


def connectivity(net, measure='alpha'):
    """Connectivity Indices
    
    Levinson D. (2012) Network Structure and City Size.
                PLoS ONE 7(1): e29721. doi:10.1371/journal.pone.0029721
    
    * alpha
    The alpha index is the ratio of the actual number of circuits
    in a network to the maximum possible number of circuits on that
    planar network. Values of a range from 0 percent - no circuits -
    to 100 percent - a completely interconnected network.
        alpha = e - v + p / 2*v - 5
    
    * beta
    The beta index measures the connectivity relating the number of
    edges to the number of nodes. The greater the value of beta,
    the greater the connectivity.
        beta = e / v
    
    * gamma
    The gamma index measures the connectivity in a network. It is a
    measure of the ratio of the number of edges in a network to the
    maximum number possible in a planar network. Gamma ranges from 0
    (no connections between nodes) to 1.0 (the maximum number of
    connections, with direct links between all the nodes).
        gamma = e / 3(v-2)
    
    * eta
    The eta index measures the length of the graph over
    the number of edges.
        eta = L(G) / e
    
    e = number of edges (links)
    v = number of vertices (nodes)
    p = number of graphs or subgraphs, and for a network where every
        place is connected to the network p = 1
    L(G) = total length of the graph
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    measure : str
        statistic to calculate
    
    Returns
    -------
     : measure in question
    """
    
    e = float(net.n_segm)
    v = float(net.n_node)
    p = float(net.n_edge_cc)
    L = net.network_length
    
    if measure == 'alpha':
        return (e-v+p) / ((2*v) - 5)
    
    if measure == 'beta':
        return e / v
    
    if measure == 'gamma':
        # number of edges in a maximally connected planar network
        e_max = (3 * (v-2))
        return e / e_max
    
    if measure == 'eta':
        return L / e


def entropy(net):
    """Entropy
    
    Levinson D. (2012) Network Structure and City Size.
                PLoS ONE 7(1): e29721. doi:10.1371/journal.pone.0029721

    entropy measure of heterogeneity
    H = -sum([pi np.log2 * (pi) for pi in subset_porportions])
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    
    Returns
    -------
    indiv_entropies : dict
        indivual entropies of MTFCC categories
    """
    
    indiv_entropies = {}
    
    for mtfcc in net.s_data.MTFCC.unique():
        
        subset = net.s_data[net.s_data.MTFCC == mtfcc].shape[0]
        
        subset_proportion = float(subset) / float(net.n_segm)
        
        entropy = subset_proportion * np.log2(subset_proportion)
        
        indiv_entropies[mtfcc] = entropy
    
    return indiv_entropies


def _get_dia_rad(mtx, stat):
    """get the diameter or radius of a network with associated nodes
    
    Parameters
    ----------
    mtx : numpy.ndarray
        cost matrix
    stat : str
        min or max. Default is 'max.'
    
    Returns
    -------
     : [idx,value]
    """
    
    if stat == 'max':
        value = mtx.max()
    else:
        value = mtx[mtx != 0.].min()
    
    idx = tuple(np.where(mtx==value))
    
    if len(idx[0]) > 1:
        idx = tuple(idx[0])
    else:
        idx = tuple([idx[0][0], idx[1][0]])
    
    return [idx,value]


def get_roots(adj):
    """create a rooted object that stores connected components
    https://stackoverflow.com/questions/10301000/
    python-connected-components
    
    Parameters
    ----------
    adj : list
        record of adjacency
    
    Returns
    -------
    ccs : list
        rooted connected components
    """
    
    
    def _find_root_depth(obj, root):
        """find the root an
        
        Parameters
        ----------
        obj : int
            index of object
        root : list
            root lookup
        
        Returns
        -------
        obj : int
            root of the object
        root[obj][1] : int
            depth of the rooted object
        """
        
        while obj != root[obj][0]:
            obj = root[obj][0]
        return obj, root[obj][1]
    
    if type(adj) == dict:
        adj = [[idx, list(cc)] for idx, cc in list(adj.items())]
    
    # 1. set all objects within the root lookup to zero
    root = {i: (i, 0) for (i, neighs) in adj}
    
    # 2. iterate through each combination of neighbors
    for (i, neighs) in adj:
        for j in neighs:
            
            # 2-A. find the root of i and its depth
            root_of_i, depth_of_i = _find_root_depth(i, root)
            
            # 2-B. find the root of j and its depth
            root_of_j, depth_of_j = _find_root_depth(j, root)
            
            # 2-C. set each object as either root or non root
            if root_of_i != root_of_j:
                _min = root_of_i
                _max = root_of_j 
                if  depth_of_i > depth_of_j:
                    _min = root_of_j
                    _max = root_of_i
                root[_max] = _max, max(root[_min][1]+1,
                                       root[_max][1])
                root[_min] = (root[_max][0],-1)
    
    # 3. initialze connected components by root lookup
    ccs = {}
    
    # 4. create empty list entry for each rooted connected compponent
    for (i, neighs) in adj:
        if root[i][0] == i:
            ccs[i] = []
    
    # 5. fill each list with the components
    [ccs[_find_root_depth(i, root)[0]].append(i) for (i, neighs) in adj]
    ccs = [list(cc) for cc in list(ccs.items())]
    
    return ccs


def get_neighbor_distances(net, v):
    """create a dictionary of neighbors with distances for dijkstra
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    v  : int
        source node for iteration of neighbors
    
    Returns
    -------
    neighbors : dict
        dictionary of {neighor ID: neighbor distance}
    """
    
    neighbors = {}
    
    for e in net.node2segm[v][1]:
        
        if net.segm2node[e][1][0] != v:
            neighbors[net.segm2node[e][1][0]] = net.segm2len[e][1]
        
        else:
            neighbors[net.segm2node[e][1][1]] = net.segm2len[e][1]
    
    return neighbors


def generate_tree(pred):
    """Generate a tree for shortest path between
    source and destination nodes
    
    Parameters
    ----------
    pred : list
        list of predecessor nodes
   
   Returns
    -------
    tree : dict
        shortest path tree in the form -- 
        {source node: [destination node ... source node]}
        e.g.: {0: [0], 1: [0], 2: [339, 0]}
    """
    
    tree = {}
    
    for i, p in enumerate(pred):
        if p == -1:
            #root node
            tree[i] = [i]
            continue
        
        idx = p
        path = [idx]
        
        while idx >= 0:
            nextnode = pred[idx]
            idx = nextnode
            if idx >= 0:
                path.append(nextnode)
        tree[i] = path
    
    return tree


def _euc_dist(p1, p2):
    """Calculate the euclidean distance between two line endpoints
    
    Parameters
    ----------
    p1 : float or int
        enpoint 1 of a line
    p2 : float or int
        enpoint 2 of a line
    
    Returns
    -------
    euclidean distance between two line endpoints
    """
    
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def shortest_path(net, gp=False):
    """graph traversal for shortest path
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
        generate paths. Default is False.
    gp : bool
        generate paths. Default is False.
    
    Returns
    -------
    mtx : numpy.ndarray
        shortest path costs between all nodes
    paths : dict
        graph traveral paths. ONLY FOR DIJKSTRA CURRENTLY.
    """
    
    # Instantiate empty cost matrix and paths
    mtx, paths = np.empty((net.n_node, net.n_node)), {}
    
    # Dijkstra
    if net.n2n_algo == 'dijkstra':
        
        # Classic source-to-all algo for optimal
        # shortest path graph traversal.
        for n in net.n_ids:
            
            # get the distance array and predecessor nodes for each node
            dist, pred = dijkstra(net, n)
            tree = None
            
            # if recording the paths
            if gp:
                tree = generate_tree(pred)
            
            # set the distance array in a matrix and paths in a dict
            mtx[n], paths[n] = dist, tree
    
    return mtx, paths


def dijkstra(net, source):
    """Dijkstra single source to all destinations
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    source : int
        source node for iteration
    
    Returns
    -------
    distance : list
        list of distances from the source node
    pred : list
        list of predecessor nodes. 
    """
    
    initial_dist = np.inf
    distance = [initial_dist for n in net.n_ids]
    distance[source] = 0.
    unvisited, pred = set([source]), [-1 for n in net.n_ids]
    
    while unvisited:
        
        # Get node with the lowest value from distance.
        dist = initial_dist
        for node in unvisited:
            if distance[node] < dist:
                dist, current = distance[node], node
        
        # Remove that node from the set.
        unvisited.remove(current)
        
        # Get the neighbors & distances to the current node.
        neighbors = get_neighbor_distances(net, current)
        for neigh, add_dist in list(neighbors.items()):
            
            new_dist = distance[current] + add_dist
            
            if distance[neigh] > new_dist:
                distance[neigh] = new_dist
                pred[neigh] = current
                unvisited.add(neigh)
    
    return distance, pred


def _check_symmetric(a, tol=1e-8):
    """validate matrix symmetry for nXn matrices
    
    Parameters
    ----------
    a : numpy.ndarray
        cost matrix
    tol : float
        tolerance. Default is 1e-8.
    
    Returns
    -------
     : bool
        symmetric or not
    """
    
    return np.allclose(a, a.T, atol=tol)


def _locate_naps(net):
    """Locate all non articulation points in order to simplfy graph
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    
    Returns
    -------
    napts : dict
        dictionary of non-articulation points and segments
    """
    
    # subset only degree-2 nodes
    degree_two_nodes = set([n for (n,d) in net.node2degree if 2 in d])
    
    # recreate n2n xwalk
    new_n2n = {k:v for (k,v) in net.node2node}
    two2two = {k: new_n2n[k] for k in degree_two_nodes}
    
    # get set intersection of degree-2 node neighbors
    for k,vs in list(two2two.items()):
        two2two[k] = list(degree_two_nodes.intersection(set(vs)))
    
    # convert back to list
    two2two = [[k,vs] for k,vs in list(two2two.items())]
    
    # created rooted non-articulation nodes object
    rooted_napts = get_roots(two2two)
    napts = {}
    napts_count = 0
    for (k,v) in rooted_napts:
        napts_count += 1
        napts[napts_count] = {net.nid_name: v}
    
    # add segment info to rooted non-articulation point object
    for napt_count, napt_info in list(napts.items()):
        napt = []
        for napt_node in napt_info[net.nid_name]:
            napt.extend([i[1] for i in net.node2segm if i[0] == napt_node])
        
        # if more than one pair of segments in napt
        napt = set([seg for segs in napt for seg in segs])
        napts[napt_count].update({net.sid_name:napt})
    
    return napts


def _simplifysegs(net, na_objs):
    """drop nodes and weld bridge segments
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    na_objs : dict
        non articulation point information
    
    Returns
    -------
    net.s_data : geopandas.GeoDataFrame
        simplified segments dataframe
    """
    
    phase = 'Simplify'
    
     # for each bridge
    for na_objs_sidx, na_objs_info in list(na_objs.items()):
        
        # get the dominant SegIDX and dataframe index
        inherit_attrs_from, idx = _get_hacky_index(net, na_objs_info)
        
        # set total length to 0 and instantiate an empty segments list
        total_length = 0.
        geoms = []
        
        # add the length of each segment to total_length and add
        # the segment to geoms
        for segm in na_objs_info[net.sid_name]:
            total_length += net.s_data.loc[(net.s_data[net.sid_name] == segm),
                                            net.len_col].squeeze()
            geom = net.s_data.loc[(net.s_data[net.sid_name] == segm),
                                   net.geo_col].squeeze()
            geoms.append(geom)
        
        # take the dominant line segment id out of the `remove` list
        na_objs_info[net.sid_name].remove(inherit_attrs_from)
        
        # add new total length cell value
        net.s_data.loc[idx, net.len_col] = total_length
        
        # add new welded line segment of dominant and non-dominant lines
        welded_line = _weld_MultiLineString(geoms)
        net.s_data.loc[idx, net.geo_col] = welded_line
        
        # remove all non-dominant line segments from the dataframe
        net.s_data = net.s_data[~net.s_data[net.sid_name]\
                                    .isin(na_objs_info[net.sid_name])]
        
        # record phase
        net.s_data.loc[idx, 'Phase'] = phase
    
    return net.s_data


def _get_hacky_index(net, ni):
    """
    VERY hacky function to get back dataframe index due to 
    trouble with using -->
    df.loc[(df[sidx] == segidx), 'geometry'\
                                  = _weldMultiLineString(geoms)
    *** Not sure if bug, but seems to be interoperability issues with ...
        shapely objects and pandas-style dataframe indexing
    *** see also:
            sauce._weld_MultiLineString()
            spaghetti.SpaghettiPointPattern.snap_to_nearest\
                                           ._record_snapped_points\
                                           ._casc2point
    
    Parameters
    ----------
    net : spaghetti.SpaghettiNetwork
    ni : dict
        non articulation point information
    
    Returns
    -------
    inherit_attrs_from : int
        segment id
    idx : int
        dataframe index value
    """
    
    df, lc, sid = net.s_data, net.len_col, net.sid_name
    
    # get maximum length
    max_len = max([df.loc[(df[sid]==segm), lc].squeeze() for segm in ni[sid]])
    
    # inherit attributes from the longest segment (SegIDX)
    inherit_attrs_from = df.loc[(df[lc] == max_len), sid].squeeze()
    
    # get the df index of 'SegIDX == inherit_attrs_from' and maxLen
    dominant = (df[sid] == inherit_attrs_from)
    longest = (df[lc] == max_len)
    idx = df.loc[dominant & longest].index[0]
    
    return inherit_attrs_from, idx


def _write_paths_csv(data, filepath):
    """Write out the paths form the generated tree dictionary
    Called from within n2nCosts()
    
    Parameters
    ----------
    data : dict
        shortest path trees
    filepath : str
        path to save
    """
    
    import csv
    
    with open(filepath+'.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for k,v in list(data.items()):
            writer.writerow([k,v])


###############################################################################
############## spaghetti.SpaghettiPointPattern functionality ##################
###############################################################################


def obs2obs_costs(orig, dest=None, dist_col=None, symmetric=False,
                  network_matrix=None, from_nodes=False, wsnap_dist=None,
                  assoc_col=None, distance_type=None, xyid=None,
                  numeric_cols=None):
    """internal function to calculate a cost matrix from (n) observations
    to (m) observations
    
    Parameters
    ----------
    orig : geopandas.GeoDataFrame
        origin observations
    dest : geopandas.GeoDataFrame
        destination observations. Default is None.
    dist_col : str
        column name for distance. ONLY used when calculating from
        endpoints to endpoints. Default is None.
    from_nodes : bool
        calcualte cost matrix from network nodes only. Default is False.
    network_matrix : numpy.ndarray
        nXn network nodes cost matrix. Default is None.
    symmetric : bool
        calculate an observation nxn cost matrix. Default is False.
    assoc_col : str
        column name for network geometry snapping. Default is None.
    wsnap_dist : str
        column name to use for distance to observation from the network.
        Default is None.
    distance_type : str
        type of distance cost matrix. Default is 'network_centroid'.
        Options are 'network_pp2n' and 'euclidean'. Default is None.
    xyid : str
        string xyID column name. Default is None.
    numeric_cols : list
        columns to preprocess to ensure numeric values. Default is None.
    
    Returns
    -------
    n2m_matrix : numpy.ndarray
        nXm cost matrix
    """
    
    
    def _ensure_numeric(o, d, cols):
        """make sure dataframe columns are in proper numeric format
        
        Parameters
        ----------
        o : geopandas.GeoDataFrame
            origin observations
        d : geopandas.GeoDataFrame
            destination observations
        cols : list
            columns to preprocess to ensure numeric values
        
        Returns
        -------
        o : geopandas.GeoDataFrame
            origin observations
        d : geopandas.GeoDataFrame
            destination observations
        """
        
        # having trouble keeping the integers and floats in numeric form
        types = [int, float, np.int64, np.float64]
        
        for col in cols:
            
            if type(o[col][0]) not in types:
                o[col] = o[col].apply(lambda x: literal_eval(x))
            
            if type(d[col][0])  not in types:
                d[col] = d[col].apply(lambda x: literal_eval(x))
        
        return o, d
    
    
    def _dist_calc(o, i, isg, d, j, jsg, matrix, na='node_a',
                   nb='node_b', da='dist_a', db='dist_b'):
        """get the cheapest cost route from snapped pt1 to snapped pt 2
        
        Parameters
        ----------
        o : geopandas.GeoDataFrame
            origin observations
        i : int
            origin index
        isg : int
            segment associated with i
        d : geopandas.GeoDataFrame
            destination observations
        j : int
            destination index
        jsg : int
            segment associated with j
        matrix : numpy.ndarray
            all node to all node network cost matrix
        na : str
            right node label (may actually be to the visual left).
            Default is 'R_node'.
        nb : str
            left node label (may actually be to the visual right).
            Default is 'L_node'.
        da : str
            distance label to 'rn'. Default is 'R_dist'.
        db : str
            distance label to 'ln'. Default is 'L_dist'.
        
        Returns
        -------
        initial_dist : float or int
            distance from snapped point to snapped point
        """
        
        # get node indices
        ai, bi = o[na][i], o[nb][i]
        aj, bj = d[na][j], d[nb][j]
        
        # origin right and left distance
        ai_dist, bi_dist = o[da][i], o[db][i]
        
        # destination right and left distance
        aj_dist, bj_dist = d[da][j], d[db][j]
        
        # if observation nodes are snapped to the same segment
        if (ai, bi) == (aj, bj) and isg == jsg:
            
            # if the ORIGIN 'right' distance is greater than
            # (or equal to) the DESTINATION 'right' distance then
            # subtract the DESTINATION from the ORIGIN distance...
            if ai_dist >= aj_dist:
                initial_dist = ai_dist - aj_dist
            
            # ... otherwise subtract the 'right' ORIGIN from
            # the DESTINATION distance
            else:
                initial_dist = aj_dist - ai_dist
        
        # observation nodes are snapped to different segments
        else:
            # get all combinations of potential distance
            a2a = matrix[ai, aj]
            a2b = matrix[ai, bj]
            b2b = matrix[bi, bj]
            b2a = matrix[bi, aj]
            
            # create distance lookup dictionary
            lookup = {(ai, aj, 'ai','aj'): a2a,
                      (ai, bj, 'ai','bj'): a2b,
                      (bi, bj, 'bi','bj'): b2b,
                      (bi, aj, 'bj','aj'): b2a}
            
            # get minimum distances and associated nodes ids
            (n1, n2, pos_i, pos_j) = min(lookup, key=lookup.get)
            initial_dist = min(lookup.values())
            
            # evaulate the cheapest cost
            # if the lookup decides the 'right' node for ORIGIN and the 
            # 'right' distance is lower than the 'left' distance add the
            #'right' distance to the initial distance
            if pos_i == 'ai':
                initial_dist += ai_dist
            
            # otherwise add the 'left' distance to the initial distance
            else:
                initial_dist += bi_dist
            
            # if the lookup decides the 'right' node for DESTINATION and
            # the 'right' distance is lower than the 'left' distance
            # add the 'right' distance to the initial distance
            if pos_j == 'aj':
                initial_dist += aj_dist
            
            # otherwise add the left distance to the initial distance
            else:
                initial_dist += bj_dist
        
        return initial_dist
    
    
    def _return_coords(xyid):
        """convert string (list) xyid into a tuple of (x,y) coordinates
        
        Parameters
        ----------
        xyid : str
            string xy ID
        
        Returns
        -------
        coords : tuple
            (x,y) coordinates
        """
        
        # coerce the id into a plain string if in another text format
        # and do literal evaluation to tease out ID from list (if list)
        xyid = str(literal_eval(str(xyid))[0])
        
        # characters to split by and ignore
        chars = ['x', 'y', '']
        
        # return numeric x and y coordinate split at 'x' and 'y'
        coords = [float(c) for c in re.split('(x|y)', xyid) if c not in chars]
        coords = tuple(coords)
        
        return coords
    
    
    # set matrix style 
    if symmetric:
        dest = copy.deepcopy(orig)
    
    # instantiate empty matrix
    n2m_matrix = np.zeros((orig.shape[0], dest.shape[0]))
    
    # Euclidean observation nodes distance matrix
    if distance_type == 'euclidean':
        for ix in orig.index:
            for jx in dest.index:
                i, j = orig[xyid][ix], dest[xyid][jx]
                p1, p2 = _return_coords(i), _return_coords(j)
                n2m_matrix[ix,jx] = _euc_dist(p1, p2)
    
    # network-style cost matrices
    else:
        
        # make sure node indices and distances are set to
        # numeric values in dataframe columns
        orig, dest = _ensure_numeric(orig, dest, numeric_cols)
        
        # Network (from 'network nodes') distance matrix
        if from_nodes:
            
            for i in orig.index:
                
                for j in dest.index:
                    
                    # if i and j are the same observation
                    # there is no distance
                    if i == j and symmetric:
                        n2m_matrix[i,j] = 0.
                    
                    else:
                        I = orig.loc[i, assoc_col]  # i network node ID
                        J = dest.loc[j, assoc_col]  # j network node ID
                        # network i to j dist
                        net_dist = network_matrix[I,J]
                        
                        # add in distance from observation to network
                        if wsnap_dist:
                            # snapped i dist
                            from_i = orig[wsnap_dist][i]
                            # snapped j dist
                            from_j = dest[wsnap_dist][j]
                            net_dist = from_i+from_j+net_dist
                        
                        n2m_matrix[i,j] = net_dist
        
        # Network (from snapped point) distance matrix
        else:
            
            # complete distance
            # start p1 --> snap point --> nearest node -->
            # furtherst(closest) node --> snap point --> goal p2
            for i in orig.index:
                
                for j in dest.index:
                    
                    # if i and j are the same observation
                    # there is no distance
                    if i == j and symmetric:
                        n2m_matrix[i,j] = 0.
                    
                    else:
                        # network segments associated with i and j
                        isegm, jsegm = orig[assoc_col][i], dest[assoc_col][j]
                        
                        # calculate network distance from i to j
                        network_dist = _dist_calc(orig, i, isegm,
                                                  dest, j, jsegm,
                                                  network_matrix)
                        
                        # add in distance from observation to network
                        if wsnap_dist:
                            orig_snap = orig[wsnap_dist][i]
                            dest_snap = dest[wsnap_dist][j]
                            snap_dist = orig_snap + dest_snap
                            network_dist += snap_dist
                        
                        n2m_matrix[i,j] = network_dist
    
    return n2m_matrix


def remove_restricted(df, net, restr=None, col=None):
    """subset the unrestricted segments and update adjacencies and
    return (1) the geodataframe of unrestricted segments, and (2) the
    network object with updated adjacency values.
    
    Used for snapping points to appropriate network segments
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        network segments dataframe. Default is None.
    net : spaghetti.SpaghettiNetwork
    restr : list
        restricted segment types. Default is None.
    col : str
        column name for segment restriction stipulation.
        Default is None.
    
    Returns
    -------
    unrestrict_segm__df : geopandas.GeoDataFrame
        network segments of unrestricted edges dataframe
    net : spaghetti.SpaghettiNetwork
    """
    
    
    def _convert_to_dict(x2x):
        """convert the list representations of adjacency to dictionary
        representations. This is done to preserve the 'id' when removing
        segments/nodes before building the KDTree. By default the KDTree
        gives a [0,1,2, ... n] id which creates a problem when
        considering the removal of segment/nodes due to the adjacency
        relationships and KDTree not having the same 'n'.
        
        Parameters
        ----------
        x2x : list
            list representation of (element)2(element) adjacency
        
        Returns
        -------
        converted_dict : dict
            dict representation of (element)2(element) adjacency
        """
        
        converted_dict = {}
        for (k,v) in x2x:
            converted_dict[k] = v
        
        return converted_dict
    
    
    def _remove_segm_ids(n2s, restr):
        """remove restricted segment ids from node to segment adjacency
        
        Parameters
        ----------
        n2s : dict
            node to segment adjacency
        restr : list
            restricted segment ids
        
        Returns
        -------
        n2s : dict
            updated node to segment adjacency
        """
        
        update_n2s = {}
        
        for node, segms in list(n2s.items()):
            updated_segms = []
            
            for segm in segms:
                if segm not in restr:
                    updated_segms.append(segm)
            update_n2s[node] = updated_segms
        
        n2s = update_n2s
        
        return n2s
    
    
    # records all node and segment ids
    full_node_ids = net.n_ids
    full_segm_ids = net.s_ids
    
    # get restricted segment ids
    restricted_segm_ids = list(df[df[col].isin(restr)].index)
    
    # subset unrestricted segments dataframe
    unrestrict_segm__df = df[~df[col].isin(restr)]
    
    # remove restricted segment ids from the segment id list
    net.s_ids = list(unrestrict_segm__df.index)
    
    # update segm2node & convert to dict
    net.segm2node = remove_adj(net.segm2node, restricted_segm_ids)
    net.segm2node = _convert_to_dict(net.segm2node)
    
    # update segm2geom & convert to dict
    net.segm2geom = remove_adj(net.segm2geom, restricted_segm_ids)
    net.segm2geom = _convert_to_dict(net.segm2geom)
    
    # convert segm2len to dict
    net.segm2len = _convert_to_dict(net.segm2len)
   
    # remove restricted node ids from the node id list
    unrestrict_node_ids = set()
    for seg, nodes in list(net.segm2node.items()):
        unrestrict_node_ids.update(nodes)
    net.n_ids = list(unrestrict_node_ids)
    
    # record restricted nodes: these are nodes that fall ONLY on
    # restricted segments (e.g. node incident with 2 restricted
    # segments). Node that are incident with 1 restricted and 1
    # unrestricted segment are NOT restricted but the node2segm
    # adjacency for the restricted segment will be removed in below.
    restricted_node_ids = [node for node in full_node_ids \
                           if node not in net.n_ids]
    
    # update node2coords & convert to dict
    net.node2coords = remove_adj(net.node2coords, restricted_node_ids)
    net.node2coords = _convert_to_dict(net.node2coords)
    
    # update node2geom & convert to dict
    net.node2geom = remove_adj(net.node2geom, restricted_node_ids)
    net.node2geom = _convert_to_dict(net.node2geom)
    
    # convert node2segm to dict
    net.node2segm = _convert_to_dict(net.node2segm)
    
    # update node2segm adjacency by removing restricted segments from
    # nodes that are incident with both restricted and unrestricted
    # segments
    net.node2segm = _remove_segm_ids(net.node2segm, restricted_segm_ids)
    
    return unrestrict_segm__df, net 


def get_obs2coords(pp):
    """create an observation to coordinate x walk
    
    Parameters
    ----------
    pp : spaghetti.SpaghettiPointPattern
    
    Returns
    -------
    o2c : list
        observations to coordinates
    """
    
    o2c = [[(ix, pp.df.loc[ix, pp.df_key]),\
            (pp.df.loc[ix, pp.geo_col].x,
             pp.df.loc[ix, pp.geo_col].y)]\
                     for ix in pp.df.index]
    
    return o2c


def snap_to_nearest(pp, sframe=None, net=None, kd_tree=None):
    """Record the nearest network node then snap to either that
    endpoint or the nearest line segment.
    
    Parameters
    ----------
    pp : spaghetti.SpaghettiPointPattern
    sframe : geopandas.GeoDataFrame
        network segments dataframe. Default is None.
    net : spaghetti.SpaghettiNetwork
        network object. Default is None.
    kd_tree  : scipy.spatial.kdtree.KDTree
        coords look up kd tree. Default is None.
    
    Returns
    -------
    snp_pts_df : geopandas.GeoDataFrame
        dataframe of newly created empirical observations
        snapped to the network.
    """
    
    
    def _get_k_nearest(pp, net):
        """ record k nearest in the form -- 
            [[(idx, key), [dists, nodes], [segments]]]
        
        Parameters
        ----------
        net : spaghetti.SpaghettiNetwork
            network object
        
        Returns
        -------
        k_nearest : list
            information on the k-nearest neighbors
        """
        
        if pp.k > len(net.s_ids):
            pp.k = len(net.s_ids)
        
        k_nearest = []
        
        for (obidx,coords) in pp.obs2coords:
            
            # query kdtree and return array-like list in the form:
            # >>>    [[dist1, n1], [dist2, node2], ...[distn, noden]]
            tree = np.array(pp.kd_tree.query(coords,
                                               k=pp.k)).T.tolist()
            
            # convert id from float to int
            tree = [[dist, net.n_ids[int(idx)]] for (dist,idx) in tree]
            
            # get all segments associated with each k-nearest node
            segms = [net.node2segm[idx] for (dist,idx) in tree]
            segms = set([seg for segs in segms for seg in segs])
            
            # get node distances and nodes from tree query
            nodes = [idx for (dist,idx) in tree]
            node_dists = [dist for (dist,idx) in tree]
            
            # append to nearest information list
            k_nearest.append([obidx, [node_dists, nodes], segms])
        
        return k_nearest
    
    
    def _record_snapped_points(pp, net, kne):
        """find the nearest point along a line segment of the nearest
        network node and records pertinent information.
        
        Parameters
        ----------
        net : spaghetti.SpaghettiNetwork
            network object
        kne : list
            information on the k-nearest neighbors
        
        Returns
        -------
        snpts : dict
            information on the newly snapped empirical observations
        """
        
        
        def _unary2point(pp, net, sf, sgs, spinfo, opt):
            """find the nearest point along the lines of a unary union
            
            Parameters
            ----------
            net : spaghetti.SpaghettiNetwork
                network object
            sf : geopandas.GeoDataFrame
                network segments dataframe
            sgs : list
                unrestricted segment ids
            spinfo : dict
                snapped point information
            opt : shapely.geometry.Point
                original empirical point geometry
            
            Returns
            -------
            spinfo : dict
                updated snapped point information
            """
            
            # get line subset
            segms_sub = sf[sf[pp.sid_name].isin(sgs)]
            segms_uu = segms_sub.unary_union
            near_points = []
            
            for line in segms_uu:
                near_pt = line.interpolate(line.project(opt))
                near_points.append(near_pt)
            
            NEAREAST_POINT_DIST = np.inf
            point = None
            for near_p in near_points:
                if near_p.distance(opt) < NEAREAST_POINT_DIST:
                    NEAREAST_POINT_DIST = near_p.distance(opt)
                    point = near_p
            
            # get coords of new point
            spinfo[pp.geo_col] = point
            
            # tease out which exact segment
            for ns in sgs:
                if net.segm2geom[ns].intersects(point.buffer(pp.tol)):
                    spinfo['assoc_segm'] = ns
                    break
            return spinfo
        
        
        # snapped points dictionary
        snpts = {}
        
        for (idx,key),(dists,nodes),(segms) in kne:
            
            spinfo = {pp.df_key: key}
            if pp.snap_to == 'segments':
                
                # original point geometry
                opt = pp.df[pp.geo_col][idx]
                spinfo = _unary2point(pp, net, sframe,
                                      segms, spinfo, opt)
                
                # associated segment
                asc_seg = spinfo['assoc_segm']
                
                # snapped point geometry
                snp = spinfo[pp.geo_col]
                line_length = net.segm2len[asc_seg]
                
                # line segment geometry
                asc_seg_geom = net.segm2geom[asc_seg]
                
                # endpoint of the associated line segment
                line_ep1 = net.segm2node[asc_seg][0]
                line_ep2 = net.segm2node[asc_seg][1]
                
                # line segment endpoint 1 geometry
                line_ep1_geom = net.node2geom[line_ep1]
                
                # set line endpoints 1 and 2 as node a and node b
                # node a is determined as being 0.0 distance from the
                # projected distance along the line segment
                if asc_seg_geom.project(line_ep1_geom) == 0.0:
                    node_a = line_ep1
                    node_b = line_ep2
                else:
                    node_a = line_ep2
                    node_b = line_ep1
                
                # node a information ----------------------------------
                #       * node a is the point which shapely chooses as
                #         the base from which calculation occurs
                #         for 'project' distance. 
                #       * node a is generally the 'smaller' of the 
                #         the two line segment endpoint in
                #         euclidean space
                dist_a = asc_seg_geom.project(snp)
                spinfo['dist_a'] = dist_a
                spinfo['node_a'] = node_a
                # node b information ----------------------------------
                dist_b = line_length - dist_a
                spinfo['dist_b'] = dist_b
                spinfo['node_b'] = node_b
                spinfo['dist2line'] = snp.distance(opt)
            
            # just to the nearest network vertex
            # does not currently stipulate tha ti has to be the
            # nearest vertex on the nearest line
            # can add in functionality later
            if pp.snap_to == 'nodes':
                nearest_node, nearest_dist = nodes[0], dists[0]
                spinfo['assoc_node'] = nearest_node
                spinfo['dist2node'] = nearest_dist
                spinfo[geo_col] = net.node2geom[nearest_node]
            
            # add to dictionary
            snpts[idx] = spinfo
        
        return snpts
    
    # within `remove_restricted()` this is converted to a dictionary
    # so convert back to dict if not been already
    if type(net.node2segm) == list:
        net.node2segm = _convert_to_dict(net.node2segm)
    if type(net.segm2geom) == list:
        net.segm2geom = _convert_to_dict(net.segm2geom)
    
    k_near_elems = _get_k_nearest(pp, net)
    
    # snap points
    snapped_pts = _record_snapped_points(pp, net, k_near_elems)
    
    # make columns headers based on snap method
    if pp.snap_to == 'segments':
        cols = [pp.df_key, 'assoc_segm', 'dist2line',
                'dist_a', 'node_a', 'dist_b','node_b', pp.geo_col]
    
    if pp.snap_to == 'nodes':
        cols = [pp.df_key, 'assoc_node', 'dist2node', pp.geo_col]
    
    # add decision variables if not incident calls
    if pp.df_name != 'ResidentIncidents'\
    and pp.study_area != 'Test_Grid_Leon_FL':
        try:
            dv = 'desc_var'
            cols = cols + [dv]
            [info_dict.update({dv: pp.df.loc[idx, dv]}) \
                                        for idx, info_dict in snapped_pts.items()]
        except KeyError:
            pass
    
    # create dataframe
    snp_pts_df = utils.fill_frame(pp.df, full=True,
                                  col=cols, data=snapped_pts)
    
    # add xyid
    node2xyid = utils.generate_xyid(df=snp_pts_df, geom_type='node',
                                    geo_col=pp.geo_col)
    snp_pts_df = utils.fill_frame(snp_pts_df, idx='index',
                                  col=pp.xyid, data=node2xyid)
    
    # create observation - to - segment lookup
    pp.obs2segm = dict(zip(snp_pts_df[pp.df_key],
                             snp_pts_df['assoc_segm']))
    
    return snp_pts_df


def _get_lines(hspace, vspace, withbox, bounds, hori=True):
    """generate line segments for a grid

    Parameters
    ----------
    hspace : list
        horizontal spacing
    vspace : list
        vertical spacing
    withbox : bool
        include outer rim
    bounds : list
        area bounds in the form of [x1,y1,x2,y2].
    hori : bool
        generate horizontal line segments.
        Default is True. False generates vertical segments.  
    
    Returns
    -------
    lines : list
        All vertical or horizontal line segments in the grid.
    """

    # Initialize starting and ending horizontal indices
    h_start_at, h_end_at = 0, len(hspace)

    # Initialize starting and ending vertical indices
    v_start_at, v_end_at = 0, len(vspace)

    # set inital index track back to 0
    y_minus = 0
    x_minus = 0

    if hori: # start track back at 1 for horizontal lines
        x_minus = 1
        if not withbox: # do not include borders
            v_start_at += 1
            v_end_at -= 1

    else: # start track back at 1 for vertical lines
        y_minus = 1
        if not withbox: # do not include borders
            h_start_at += 1
            h_end_at -= 1

    # Create empty line list and fill
    lines = []

    # for element in the horizontal index
    for hplus in range(h_start_at, h_end_at):

        # for element in the vertical index
        for vplus in range(v_start_at, v_end_at):

            # ignore if a -1 index
            if hplus-x_minus == -1\
            or vplus-y_minus == -1:
                continue
            else:
                # Point 1 (start point + previous slot in 
                #          horizontal or vertical space index)
                p1x = bounds[0] + hspace[hplus-x_minus]
                p1y = bounds[1] + vspace[vplus-y_minus]
                p1 = Point(p1x, p1y)

                # Point 2 (start point + current slot in 
                #          horizontal or vertical space index)
                p2x = bounds[0] + hspace[hplus]
                p2y = bounds[1] + vspace[vplus]
                p2 = Point(p2x, p2y)

                # LineString
                lines.append(LineString((p1,p2)))
    return lines


def _polygonize(geoms):
    """polygonize shapely.LineString objects
    
    Parameters
    ----------
    geoms : list
        contiguous LineString objects forming a LinearRing
    
    Returns
    -------
    geoms : shapely.Polygon
        polygonized LineString objects
    """
    geoms = list(polygonize(geoms))
    return geoms


