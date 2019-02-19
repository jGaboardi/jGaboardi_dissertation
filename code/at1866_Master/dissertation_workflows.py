"""
Disseration research Workflows for PP2N cost matrix preparation

    - for snapping points to networks
    - for calculating cost mstrices
    - performing euclidean method
    - performing the line voronoi method
    - performing the PP2N method

"""

# standard library imports
import os, time

# non-standard library imports
import geopandas as gpd
import numpy as np
import pulp

# project library imports
from . import utils
from . import spaghetti as spgh
from . import allocation as alct


def generate_grid(area, place_time, initial, c_cen, as_fixed=True,
                  bounds=None, hl=None, vl=None, fixed_population=None,
                  xyid=None, geo_col=None, sid_name=None, proj_init=None,
                  mtfcc='MTFCC', mtfcc_label='S1400', edge_file_prefix='Edges',
                  synth_id='GEOID_syn', synth_pop='synth_pop',
                  synth_poly_title='CensusBlocks',
                  synth_polypop_title='Populated_CensusBlocks',
                  synth_centroid_title='Populated_Centroid_CensusBlocks',
                  synth_hh='households', file_type='.shp',
                  synth_hh_title='Households_Synthetic'):
    '''Generate and write out synthetic grid-based data
    
    Parameters
    ----------
    area : str
        see `STUDY_AREA` in `runner_dissertation.py`
    place_time : str
        see `place_time` in `runner_dissertation.py`
    initial : str
        see `initial` in `runner_dissertation.py`
    c_cen : str
        see `c_cen` in `runner_dissertation.py`
    as_fixed : bool
        when [True] use the simple test fir disseration.
    bounds : list
        grid lines bounding box in the form of [x1,y1,x2,y2].
    hl : int
        Count of horizontal lines. Default is None.
    vl : int
        Count of vertical lines. Default is None.
    fixed_population : list
        population count per synthetic block
    xyid : str
        string xyID column name. Default is None.
    geo_col : str
        geometry column name. Default is None.
    sid_name : str
        segment column name. Default is None.
    proj_init : int
        intial coordinate reference system. default is None.
    mtfcc : str
        MTFCC dataframe column name. Default is 'MTFCC'.
    mtfcc_label : str
        feature class code. Default is 'S1400'.
    edge_file_prefix : str
        write out edge .shp file prefix. Default is 'Edges'.
    synth_id : str
        synthetic data ids. Default is 'GEOID_syn'.
    synth_pop : str
        synthetic data weights. Default is 'synth_pop'.
    synth_poly_title : str
        synthetic polygons .shp file name.
        Default is 'CensusBlocks'.
    synth_polypop_title : str
        synthetic populated polygons .shp file name.
        Default is 'Populated_CensusBlocks'.
    synth_centroid_title : str
        synthetic ppolygon centroids .shp file name.
        Default is 'Populated_Centroid_CensusBlocks'.
    synth_hh : str
        synthetic household data points. Default is 'households'.
    synth_hh_title : str
        synthetic households .shp file name.
        Default is 'Households_Synthetic'.
    file_type : str
        file extension. Default is '.shp'.
    '''
    
    if as_fixed:
        # set grid parameters
        bounds = [0,0,9,9]
        hl, vl = 2, 2
        fixed_population = [0,0,3,3,6,6,9,9,12]
    
    # create basic grid arcs
    grid_arcs = spgh.synthetic_grid(bounds=bounds,
                                    n_hori_lines=hl, n_vert_lines=vl,
                                    withbox=False, as_polys=False)
    grid_arcs = utils.set_crs(grid_arcs, proj_init=proj_init)
    grid_arcs[sid_name] = grid_arcs.index
    grid_arcs[mtfcc] = [mtfcc_label] * grid_arcs.shape[0]
    edges = edge_file_prefix + place_time
    grid_edge_dir = initial + edges
    if not os.path.exists(grid_edge_dir):
        os.makedirs(grid_edge_dir)
    grid_arcs.to_file(grid_edge_dir + '/' + edges + file_type )
    
    # create basic grid polygons
    grid_polys = spgh.synthetic_grid(bounds=bounds,
                                     n_hori_lines=2, n_vert_lines=2,
                                     withbox=True, as_polys=True)
    grid_polys = utils.set_crs(grid_polys, proj_init=proj_init)
    grid_polys[synth_id] = grid_polys.index.astype(int)
    grid_polys[synth_pop] = fixed_population
    # all polygons
    polys = synth_poly_title + place_time
    grid_polys.to_file(c_cen + polys + file_type)
    # populated only
    grid_polys = grid_polys[grid_polys[synth_pop] > 0]
    polys = synth_polypop_title + place_time
    grid_polys.to_file(c_cen + polys + file_type)
    
    # create basic grid polygon centroids
    grid_polys_centroids = grid_polys.copy()
    grid_polys_centroids.geometry = grid_polys.centroid
    # add xyid
    cent2xyid = utils.generate_xyid(df=grid_polys_centroids, geo_col=geo_col)
    idx_adjust = len(fixed_population) - grid_polys_centroids.shape[0]
    grid_polys_centroids = utils.fill_frame(grid_polys_centroids,
                                            col=xyid, data=cent2xyid,
                                            add_factor=idx_adjust)
    cents = synth_centroid_title + place_time
    grid_polys_centroids.to_file(c_cen + cents + file_type)
    
    # create synthetic households
    grid_hh = utils.create_synthetic_locs(None, area=area, data_type=synth_hh)
    grid_hh = utils.set_crs(grid_hh, proj_init=proj_init)
    if as_fixed:
        grid_hh[synth_id] = [8, 7, 6, 5, 4, 3, 2]
    
        # align population
        grid_hh[synth_pop] = [0] * grid_hh.shape[0]
        for poly in grid_polys[synth_id]:
            if poly in list(grid_hh[synth_id]):
                pop = grid_polys.loc[(grid_polys[synth_id] == poly),\
                                                 synth_pop].squeeze()
                grid_hh.loc[(grid_hh[synth_id] == poly), synth_pop] = pop
    
    # add xyid
    hh2xyid = utils.generate_xyid(df=grid_hh, geo_col=geo_col)
    grid_hh = utils.fill_frame(grid_hh, col=xyid, data=hh2xyid)
    syth_hh = synth_hh_title + place_time
    grid_hh.to_file(c_cen + syth_hh + file_type )


def allocate(network, allocate_methods, segm_file=None, node_file=None,
             clean=None, net_dir=None, cen_dir=None, inter=None,
             alloc_dir=None, geographic_units=None, restrict=None,
             geo_col=None, sid_name=None, xyid=None, desc_var=None,
             proj_init=None, pp2n_offset=10., vor_rho=300, vor_offset=2.,
             restrict_col='SegMTFCC', mtfcc='MTFCC', file_type='.shp'):
    """top-level function caller for allocation population to networks.
    Creates and writes out GeoDataFrame objects for `pp2n` and `va2n` methods.
    
    Parameters
    ----------
    network : spaghetti.SpaghettiNetwork
    allocate_methods : list
        methods for population to network allocation.
    segm_file : str
        path and file name for segment data. Default is None.
    node_file : str
        path and file name for segment data. Default is None.
    clean : str
        path to `clean` data. Default is None.
    net_dir : str
        path to `network` data. Default is None.
    cen_dir : str
        path to `census` data. Default is None.
    inter : str
        file path to intermediary data. Default is None.
    alloc_dir : str
        path to `allocation` data. Default is None.
    geographic_units : list
        census geography areas to run. Default is None.
    restrict : list
        no population on roads that should never have population.
        see `SNAP_RESTRICT` in `runner_dissertation.py`.
    geo_col : str
        geometry column name. Default is None.
    sid_name : str
        segment column name. Default is None.
    desc_var : str
        decision variable column name. Default is None.
    xyid : str
        combined x-coord + y-coords string ID. Default is None.
    proj_init : int
        intial coordinate reference system. default is None.
    pp2n_offset : {float, int}
        pp2n parallel outset (meters). Default is 10.
    vor_offset : {float, int}
        'voronoi offset'. distance from articulation point to begin
        and end segment densification. Default is 2. (meters).
    vor_rho : int
        'voronoi rho'. intial count of points to generate along line
        segments. Default is 300.
    restrict_col : str
        segment MTFCC label. Default is 'SegMTFCC'.
    mtfcc : str
        MAF/TIGER Feature Class code. Default is 'mtfcc'.
    file_type : str
        file extension. Default is '.shp'.
    """
    
    if not segm_file:
        segm_file = net_dir + 'SimplifiedSegms' + network.place_time
    net_segms = gpd.read_file(segm_file+file_type)
    net_segms = utils.set_crs(net_segms, proj_init=proj_init)
    
    if not node_file:
        node_file = net_dir + 'SimplifiedNodes' + network.place_time
    net_nodes = gpd.read_file(node_file+file_type)
    net_nodes = utils.set_crs(net_nodes, proj_init=proj_init)
    
    # get observation types dictionary
    observation_types = obsv_types(clean=clean, place_time=network.place_time)
    observation_types = observation_types['network']
    
    # for `pp2n` and LVD/Overlay
    for method in allocate_methods:
        
        # instantiate dictionary to house created gdfs and outfile names
        gdf_label = 'gdf'
        out_file_label = 'outfile'
        add_dv_write = {}
        # instantiate xyid set container
        # -- see `location_based_dv` for explanation
        xyid_set = set()
        
        for gu in geographic_units:
            
            # Start timer
            time_start = time.time()
            
            # geography unit information
            geo_unit_info = observation_types[method][gu]
            
            # path and filename to polygon .shp 
            path = geo_unit_info['poly_path']
            file_ = geo_unit_info['poly_file']
            path_file_ = path + file_ + file_type
            
            # path and filename to write out .shp
            if method == 'pp2n':
                out_file_name = '%s_%s_%s' % (method, pp2n_offset, file_)
            elif method == 'va2n':
                out_file_name = '%s_%s_%s' % (method, vor_rho, file_)
            
            out_file = alloc_dir + out_file_name + file_type
            
            # if the file already exists skip
            out_file_exists = os.path.exists(out_file)
            if not out_file_exists:
                
                # polygon file to base pp2n or va2n on
                gdf = gpd.read_file(path_file_)
                gdf = utils.set_crs(gdf, proj_init=proj_init)
                
                # polygon key and population column names
                poly_key = geo_unit_info['poly_key']
                poly_pop = geo_unit_info['poly_pop']
                
                if network.study_area == 'Test_Grid_Leon_FL':
                    poly_key = poly_key + '_syn'
                    poly_pop = 'synth_pop'
                
                if network.study_area != 'Test_Grid_Leon_FL':
                    # drop these columns before spatial join
                    geog_columns = ['AREALAND', 'AREAWATER',
                                    'CENTLAT', 'CENTLON', 'FUNCSTAT',
                                    'INTPTLAT', 'INTPTLON',
                                    'LSADC', 'OBJECTID', 'OID',
                                    'STGEOMETRY', 'STGEOMET_1',
                                    'UR', desc_var, xyid]
                    try:
                        gdf['LWBLKTYP']
                        geog_columns.append('LWBLKTYP')
                    except KeyError:
                        pass
                    gdf.drop(geog_columns, axis=1, inplace=True)
                
                # --- generate pp2n points ---
                if method == 'pp2n':
                    
                    # columns for produced dataframe
                    columns = [geo_col, restrict_col, sid_name]
                    
                    drop_columns = ['index_right']
                    
                    # call pp2n
                    gdf = alct.pp2n(net_segms, gdf, columns=columns,
                                    offset_=pp2n_offset, remove_segm=restrict,
                                    drop_cols=drop_columns, pop_type=float,
                                    sid_name=sid_name, geo_col=geo_col,
                                    xyid=xyid, poly_key=poly_key,
                                    poly_pop=poly_pop, proj_init=proj_init,
                                    restrict_col=restrict_col)
                
                # --- generate va2n points ---
                if method == 'va2n':
                    
                    # get envelope bounds and bounds to clip by
                    if network.study_area == 'Test_Grid_Leon_FL':
                        bound_type = 'CensusBlocks' + network.place_time
                        bounds_file = cen_dir + bound_type + file_type
                    else:
                        bounds_info = observation_types[method]['tracts']
                        bound_path = bounds_info['poly_path']
                        bound_file_ = bounds_info['poly_file']
                        bounds_file = bound_path + bound_file_ + file_type
                    
                    # all geographies .shp
                    _bounds = gpd.read_file(bounds_file)
                    # unary union for clipping
                    clip_by = gpd.GeoDataFrame(geometry=[_bounds.unary_union])
                    clip_by = utils.set_crs(clip_by, proj_init=proj_init)
                    # bounding envelope for initial voronoi trim
                    clip_by_env = [clip_by.envelope.squeeze()]
                    bounds = gpd.GeoDataFrame(geometry=clip_by_env)
                    bounds = utils.set_crs(bounds, proj_init=proj_init)
                    
                    if network.study_area != 'Test_Grid_Leon_FL':
                        cols_in = ['STATEFP', 'COUNTYFP', 'TLID',
                                    sid_name, mtfcc]
                    else:
                        cols_in = [sid_name, mtfcc]
                    
                    # derivation from census geography file name
                    va2n_polys = geo_unit_info['file']
                    
                    # call va2n
                    gdf = alct.va2n(network.study_area, alloc_dir, gdf,
                                    net_segms, net_nodes, inter=inter,
                                    sid_name=sid_name, clip_by=clip_by,
                                    remove_segm=restrict, xyid=xyid,
                                    bounds=bounds, cols_in=cols_in,
                                    proj_init=proj_init, poly_key=poly_key,
                                    poly_pop=poly_pop, initial_rho=vor_rho,
                                    offset_param=vor_offset,
                                    mtfcc=mtfcc, geo_col=geo_col,
                                    va2n_polys=va2n_polys)
                
                # update gdf dict
                add_dv_write[gu] = {gdf_label: gdf,
                                    out_file_label: out_file}
                # update high-precision xyid
                xyid_set.update(list(gdf[xyid]))
                
            time_end = round((time.time()-time_start)/60., 3)
            print('\t\tAllocation - %s complete in %s min.' % (method, 
                                                               time_end))
        
        # add location-based decision variable and write GeoDatFrame
        add_location_based_dv_and_write(add_dv_write, xyid_set,
                                        desc_var, xyid,
                                        gdf_label, out_file_label)


def add_location_based_dv_and_write(add_dv_write, xyid_set,
                                    desc_var, xyid,
                                    gdf_label, out_file_label,
                                    dv_type='client'):
    """create label decision variables based on the observation's
    `xyid`. This is performed in order to keep decision variable
    indexing consistent across spatial aggregation level for the
    same method of allocation (and only the same method of
    allocation). Then write the geopandas.GeoDataFrame to file.
    
    Parameters
    ----------
    add_dv_write : dict
        information lookup storing the geopandas.GeoDataFrame
        and file name for writing out.
    xyid_set : set
        all xyIDs present across spatial extents from a single
        method of allocation.
    desc_var : see `allocate()`
    xyid : see `allocate()`
    gdf_label : str
        dictionary key for the geopandas.GeoDataFrame.
    out_file_label : str
        dictionary key for the file name.
    dv_type : str
        Default to 'client' for all population-based
        allocation methods.
    """
    
    # fetch decision variables list
    dvs_list = utils.desc_var_list(dv_type, len(xyid_set))
    
    # create xyid - to - decision variable lookup
    xyid_2_dv = dict(zip(xyid_set, dvs_list))
    
    # add decision variables based on high-precision xyid
    for method, info in add_dv_write.items():
        info[gdf_label][desc_var] = info[gdf_label][xyid].apply(lambda xy:\
                                                                xyid_2_dv[xy])
    
    # write out file
    for method, info in add_dv_write.items():
        
        info[gdf_label].to_file(info[out_file_label])


def snap_obsvs(network, segm_file=None, clean=None, geographic_units=None,
               non_geographic_units=None, restrict=None, restrict_col=None,
               snap_to=None, representation=None, pp2n=None, va2n=None,
               small_unit=None, file_type='.shp'):
    """snap emprical observations to a spaghetti.Spaghetti network
    object and write out a geopandas.GeoDataFrame.
    
    Parameters
    ----------
    network : spaghetti.SpaghettiNetwork
    segm_file : str
        path to segments file
    restrict            list
        restricted segment types. Default is None.
    restrict_col : str
        column name for segment restriction stipulation.
        Default is None.
    clean : str
        file path to save clean data. Default is None.
    snap_to : str
        snap points to either segments of nodes
    representation : list
        method for allocating data to the network
    geographic_units : list
        descriptive list of geographies
    non_geographic_units :list
        descriptive list of non-census units (observations)
    pp2n : {float,int}
        offest distance for `pp2n`. Default is None.
    va2n : : {float,int}
        LVD rho for `va2n`. Default is None.
    small_unit : str
        either 'synth_households', 'households', or 'parcels'.
    file_type : str
        file extension. Default is '.shp'.
    """
    
    if not segm_file:
        segm_file = clean+'network_data/SimplifiedSegms' + network.place_time
    segmsdf = gpd.read_file(segm_file+file_type)
    
    # get observation types dictionary
    observation_types = obsv_types(clean=clean, pp2n=pp2n, va2n=va2n,
                                   place_time=network.place_time)
    observation_types = observation_types['network']
    
    if restrict:
        # segment connection restriction
        segmsdf, network = spgh.remove_restricted_segms(segmsdf, network,
                                                        restr=restrict,
                                                        col=restrict_col)
    # build kd tree
    net_nodes_kdtree = spgh.build_net_nodes_kdtree(network.node2coords)
    
    # see long note a bit below
    skip_methods = ['pp2n', 'va2n']
    
    # snap points to network for each snapping method desired
    for gu in non_geographic_units + geographic_units:
        
        if gu in non_geographic_units:
            launch_snap_job(observation_types[gu], net=network,
                            snap_to=snap_to, segm_file=segmsdf,
                            kd_tree=net_nodes_kdtree, file_type=file_type)
        
        else:
            
            for repr in representation:
                
                # These are skipped to to my set up of the
                # `observation_types` lookup. What this say doing is
                # skipping if "micro-style" units and not-centroid
                # allocation methods are employed (because the
                # "micro-style" units are not derived from either
                # the pp2n or va2n methods)
                if gu == small_unit and repr in skip_methods:
                    continue
                
                launch_snap_job(observation_types[repr][gu],
                                repr=repr, net=network,
                                snap_to=snap_to, segm_file=segmsdf,
                                kd_tree=net_nodes_kdtree,
                                file_type=file_type)


def launch_snap_job(info, repr=None, net=None, snap_to=None,
                    kd_tree=None, segm_file=None,
                    file_type=None, k=20, tol=.001):
    """launch a single observation pattern point snap process
    
    Parameters
    ----------
    info : dict
        observation points information
    repr : str
        method of representation
    kd_tree : scipy.spatial.kdtree.KDTree
        all network nodes lookup
    net : see snap_obsvs()
    snap_to : see snap_obsvs()
    k : see snap_obsvs()
    segm_file : see snap_obsvs()
    tol : see snap_obsvs()
    file_type : str
        file extension. Default is '.shp'.
    """
    
    df_name, path_to, df_key, df_file, df_pop = get_info(info)
    obsv_file = path_to + df_file + file_type
    obvs_df = gpd.read_file(obsv_file)
    
    try:
        obvs_df[df_key]
    except KeyError as error_key:
        error_key = error_key.args[0]
        if error_key == 'GEOID':
            df_key = 'GEOID_syn'
        else:
            raise KeyError
    
    repr = ''
    snapped_file = path_to + 'Snapped_' + repr + df_file + file_type
    pickle_file = 'PointPattern_' + df_name
    
    # determine existence of files
    snapped_file_exists = os.path.exists(snapped_file)
    pickle_file_exists = os.path.exists(path_to+pickle_file)
    
    if not snapped_file_exists and not pickle_file_exists:
        # instantiate pointpattern
        pointpattern = spgh.SpaghettiPointPattern(df=obvs_df, net=net,
                                                  snap_to=snap_to, k=k,
                                                  df_name=df_name, tol=tol,
                                                  df_key=df_key, df_pop=df_pop,
                                                  kd_tree=kd_tree,
                                                  net_segms=segm_file)
        # write snapped points to file
        pointpattern.snapped_points.to_file(snapped_file)
        # pickle point pattern
        spgh.dump_pickled(pointpattern, path_to, pickle_name=pickle_file)


def matrix_calc(network, segm_file=None, clean=None, snap_to=None,
                mtx_to_csv=False, df_to_csv=False, representation=None,
                geographic_units=None, non_geographic_units=None,
                pp2n=None, va2n=None, nearest_to_pickle=False,
                small_unit=None, in_rdc=False, file_type='.shp'):
    """(1) populate an origin by destination cost matrix which may be
    a numpy.ndarray written to file as .csv and/or a pandas.DataFrame
    written to file as .csv; and (2) create a nearest neighbor
    dictionary to be pickled.
    
    Parameters
    ----------
    network : spaghetti.Spaghetti
    segm_file : str
        path to segments file
    clean : str
        file path to save clean data. Default is None.
    snap_to : str
        snap points to either segments of nodes 
    df_to_csv : bool
        save the cost matrix out as a pandas.DataFrame with column and
        row observation indices. Default is False
    mtx_to_csv : bool
        save the cost matrix. Default is False.
    representation : list
        method for allocating data to the network
    pp2n : {float,int}
        offest distance for `pp2n`. Default is None.
    va2n : {float,int}
        LVD rho for `va2n`. Default is None.
    geographic_units : list
        descriptive list of geographies
    non_geographic_units : list
        descriptive list of non-census units (observations)
    nearest_to_pickle : bool
        pickle nearest neighbor dictionary (True). Default is False.
    small_unit : str
        either 'synth_households', 'households', or 'parcels'.
    in_rdc : bool
        inside the RDC [True] or not [False]. Default is False.
    file_type : str
        file extension. Default is '.shp'.
    """
    
    
    def prep_settings(snap_to):
        """internal helper function to prepare initial matrix
        calculation settings based on snap location.
        
        Parameters
        ----------
        snap_to : see parent function
        
        Returns
        -------
        from_nodes : bool
        numeric_cols : list
        wsnap : str
        assoc_col : str
        """
       
        # set snapped to nodes
        if snap_to == 'nodes':
            from_nodes = True
            numeric_cols = ['assoc_node', 'dist2node']
            wsnap = 'dist2node'
            assoc_col = 'assoc_node'
        
        # set snapped to segments
        elif snap_to == 'segments':
            from_nodes = False
            numeric_cols = ['assoc_segm', 'dist2line','dist_a',
                            'node_a', 'dist_b','node_b']
            wsnap = 'dist2line'
            assoc_col = 'assoc_segm'
        
        else:
            raise ValueError('Not able to snap to:', snap_to)
        
        return from_nodes, numeric_cols, wsnap, assoc_col
    
    
    if not segm_file:
        segm_file = clean+'network_data/SimplifiedSegms' + network.place_time
    segmsdf = gpd.read_file(segm_file+file_type)
    
    # get observation types dictionary
    observation_types = obsv_types(clean=clean, pp2n=pp2n, va2n=va2n,
                                   place_time=network.place_time)
    
    euclidean = 'euclidean'
    snapped = 'Snapped_'
    distance_types = [euclidean, 'network']
    
    # all combinations of orig - dest
    obs_combos = _get_matrix_combos(geographic_units,
                                    non_geographic_units,
                                    distance_types,
                                    representation)
    
    # see long note a bit below
    skip_methods = ['pp2n', 'va2n']
    
    for combo_info in obs_combos:
    
        if len(combo_info) == 3:
            repr, dist = combo_info[0], combo_info[1]
            orig, dest = combo_info[2], combo_info[2]
        else:
            repr, dist = combo_info[0], combo_info[1]
            orig, dest = combo_info[2], combo_info[3]
        
        # These are skipped to to my set up of the
        # `observation_types` lookup. What this say doing is
        # skipping if "micro-style" units and not-centroid
        # allocation methods are employed (because the
        # "micro-style" units are not derived from either
        # the pp2n or va2n methods)
        if dest == small_unit and repr in skip_methods:
            continue
        
        # skip if (a) census geography to census geography; and
        #         (b) county level; and
        #         (c) outside to RDC
        if orig == dest\
        and network.study_area == 'Leon_FL'\
        and not in_rdc:
            continue
        
        # set distance type, origin type,
        # destination type, and symmetric
        if orig == dest:
            symmetric = True
        else:
            symmetric = False
        
        # set shapefiles to read 'snapped' prefixes
        #  -- census geographies are snapped to the network following
        #     various methods of aggregation
        #  -- ResidentialIncidents and FireStations are simply
        #     snapped to the network as is.
        if dist != euclidean:
            _snapped = snapped
        else:
            _snapped = ''
        
        # orig observations
        if orig.startswith('F'):
                orig_info = get_info(observation_types[dist][orig])
        else:
            orig_info = get_info(observation_types[dist][repr][orig])
        obsname1, obspath1, obs_key1, obsfile1, obspop1 = orig_info
        
        # read gdf
        odsvshp1 = _snapped + obsfile1 + file_type
        odsv_df1 = gpd.read_file(obspath1+odsvshp1)
        
        # dest observations
        if symmetric:
            obsname2, obspath2, obs_key2, obsfile2, obspop2 = orig_info
        else:
            dest_info = get_info(observation_types[dist][repr][dest])
            obsname2, obspath2, obs_key2, obsfile2, obspop2 = dest_info
        
        # read gdf
        odsvshp2 = _snapped + obsfile2 + file_type
        odsv_df2 = gpd.read_file(obspath2+odsvshp2)
        
        # prepare initial matrix calculation settings based on snap
        from_nodes, numeric_cols, wsnap, assoc_col = prep_settings(snap_to)
        
        # record standard paths
        write_path = '%scost_matrices/' % clean
        name = '%s_%s_to_%s_' % (dist, obsname1, obsname2)
        write_path = '%s%s' % (write_path, name)
        
        # dataframe write path
        df_wp = write_path + 'DataFrame.csv'
        # matrix csv write path
        mtx_wp = write_path + 'CostMatrix.csv'
        # nearest neighbor pickled object write path
        nn = 'Nearest'
        nn_wp = write_path + nn
        
        # determine existence of files
        df_file_exists = os.path.exists(df_wp)
        mtx_file_exists = os.path.exists(mtx_wp)
        nn_file_exists = os.path.exists(nn_wp)
        
        if not df_file_exists\
        and not mtx_file_exists\
        and not nn_file_exists:
            
            # fill cost matrix
            mtx = spgh.obs2obs_cost_matrix(odsv_df1, dest=odsv_df2,
                                           symmetric=symmetric,
                                           network_matrix=network.n2n_matrix,
                                           from_nodes=from_nodes,
                                           wsnap_dist=wsnap,
                                           distance_type=dist,
                                           xyid=network.xyid,
                                           numeric_cols=numeric_cols,
                                           assoc_col=assoc_col)
            
            # store column and row indices to use in dataframe
            try:
                odsv_df1[obs_key1]
            except KeyError:
                obs_key1 = 'GEOID_syn'
            try:
                odsv_df2[obs_key2]
            except KeyError:
                obs_key2 = 'GEOID_syn'
            
            col_idx = [idx for idx in odsv_df1[obs_key1]]
            row_idx = [idx for idx in odsv_df2[obs_key2]]
            
            if df_to_csv:
                cix, rix = row_idx, col_idx
                frame_matrix = gpd.GeoDataFrame(mtx, columns=cix, index=rix)
                frame_matrix.to_csv(df_wp, header=True, index=True)
            
            # write out matrix with row and column indices
            if mtx_to_csv:
                np.savetxt(mtx_wp, mtx, delimiter=',')
            
            # calculate nearest neighbor step
            nearest = spgh.obs2obs_nearestneighbor(mtx, orig_ids=col_idx,
                                                   dest_ids=row_idx,
                                                   symmetric=symmetric,
                                                   keep_zero_dist=True)
            if nearest_to_pickle:
                spgh.dump_pickled(nearest, write_path, pickle_name=nn)


def facility_location():
    """
    """
    pass



def plot_geographies():
    """
    """
    pass


def plot_streets():
    """
    """
    pass


def plot_euclidean_v_network():
    """
    """
    pass


def viz_facility_location_():
    """
    """
    pass
    

def plot_facility_location():
    """
    """
    


def obsv_types(clean=None, place_time=None, pp2n='', va2n=''):
    """get dictonary of pertinent information
    regarding observation data.
    
    Parameters
    ----------
    clean : str
        path to `clean` data directory
    place_time : str
        study area within county
    pp2n : {float, int}
        label for pp2n data
    va2n : {float, int}
        label for va2n data
    
    Returns
    -------
    obs : dict
        observation information
    """
    
    euclidean_distance = 'euclidean'
    network_distance = 'network'
    
    # create euclidean entries
    obs = {euclidean_distance:
            {'pc2n':
                {'parcels':
                    {'name': 'WeightedParcels',
                     'path': clean + 'census_data/',
                     'file': 'WeightedParcels' + place_time,
                     'key': 'PARCEL_ID',
                     'pop': 'SUM_EST_PO'},
                'synth_households':
                    {'name': 'HouseholdsSynthetic',
                     'path': clean + 'census_data/',
                     'file': 'Households_Synthetic' + place_time,
                     'key': 'GEOID_syn',
                     'pop': 'POP100'},
                'pop_blocks':
                    {'name': 'pc2nPopulatedBlocks',
                     'path': clean + 'census_data/',
                     'file': 'Populated_Centroid_CensusBlocks' + place_time,
                     'key': 'GEOID',
                     'pop': 'POP100'},
                 'block_groups':
                    {'name': 'pc2nBlockGroups',
                     'path': clean + 'census_data/',
                     'file': 'Centroid_CensusBlockGroups' + place_time,
                     'key': 'GEOID',
                     'pop': 'POP100'},
                 'tracts':
                     {'name': 'pc2nTracts',
                     'path': clean + 'census_data/',
                     'file': 'Centroid_CensusTracts' + place_time,
                     'key': 'GEOID',
                     'pop': 'POP100'}},
           
            'va2n':
                {'pop_blocks':
                    {'name': 'va2nPopulatedBlocks',
                     'poly_name': 'PopulatedBlocks',
                     'path': clean + 'allocation_data/',
                     'poly_path': clean + 'census_data/',
                     'file': 'va2n_%s_Populated_CensusBlocks%s' %\
                                                                  (va2n,
                                                                   place_time),
                     'poly_file': 'Populated_CensusBlocks' + place_time,
                     'key': 'va2n_id',
                     'poly_key': 'GEOID',
                     'pop': 'pop_va2n',
                     'poly_pop': 'POP100'},
                 'block_groups':
                    {'name': 'va2nBlockGroups',
                     'poly_name': 'BlockGroups',
                     'path': clean + 'allocation_data/',
                     'poly_path': clean + 'census_data/',
                     'file': 'va2n_%s_CensusBlockGroups%s' % (va2n,
                                                               place_time),
                     'poly_file': 'CensusBlockGroups' + place_time,
                     'key': 'va2n_id',
                     'poly_key': 'GEOID',
                     'pop': 'pop_va2n',
                     'poly_pop': 'POP100'},
                 'tracts':
                     {'name': 'va2nTracts',
                      'poly_name': 'Tracts',
                      'path': clean + 'allocation_data/',
                      'poly_path': clean + 'census_data/',
                      'file': 'va2n_%s_CensusTracts%s' % (va2n, place_time),
                      'poly_file': 'CensusTracts' + place_time,
                      'key': 'va2n_id',
                      'poly_key': 'GEOID',
                      'pop': 'pop_va2n',
                      'poly_pop': 'POP100'}},
            
            'pp2n':
                {'pop_blocks':
                    {'name': 'pp2nPopulatedBlocks',
                     'poly_name': 'PopulatedBlocks',
                     'path': clean + 'allocation_data/',
                     'poly_path': clean + 'census_data/',
                     'file': 'pp2n_%s_Populated_CensusBlocks%s' % (pp2n,
                                                                   place_time),
                     'poly_file': 'Populated_CensusBlocks' + place_time,
                     'key': 'pp2n_id',
                     'poly_key': 'GEOID',
                     'pop': 'pop_pp2n',
                     'poly_pop': 'POP100'},
                 'block_groups':
                    {'name': 'pp2nBlockGroups',
                     'poly_name': 'BlockGroups',
                     'path': clean + 'allocation_data/',
                     'poly_path': clean + 'census_data/',
                     'file': 'pp2n_%s_CensusBlockGroups%s' % (pp2n,
                                                              place_time),
                     'poly_file': 'CensusBlockGroups' + place_time,
                     'key': 'pp2n_id',
                     'poly_key': 'GEOID',
                     'pop': 'pop_pp2n',
                     'poly_pop': 'POP100'},
                 'tracts':
                     {'name': 'pp2nTracts',
                      'poly_name': 'Tracts',
                      'path': clean + 'allocation_data/',
                      'poly_path': clean + 'census_data/',
                      'file': 'pp2n_%s_CensusTracts%s' % (pp2n, place_time),
                      'poly_file': 'CensusTracts' + place_time,
                      'key': 'pp2n_id',
                      'poly_key': 'GEOID',
                      'pop': 'pop_pp2n',
                     'poly_pop': 'POP100'}},
        
        'FireStations':
                    {'name': 'FireStations',
                     'path': clean + 'observation_data/',
                     'file': 'FireStations' + place_time,
                     'key': 'STATION_NU',
                     'pop': '--'},
        
        'FireStationsSynthetic':
                    {'name': 'FireStationsSynthetic',
                     'path': clean + 'observation_data/',
                     'file': 'FireStations_Synthetic' + place_time,
                     'key': 'STATION_NU',
                     'pop': '--'}
                     }
        }
    
    # copy over for network euclidean entries
    obs[network_distance] = obs[euclidean_distance]
    
    return obs


def get_info(info):
    """fetch observation .shp / data information
    
    Parameters
    ----------
    info : dict
        nested information about observations
    
    Returns
    -------
    info : tuple
        df_name, path_to, df_key, df_file, df_pop
    
    """
    
    df_name = info['name']
    path_to = info['path']
    df_key = info['key']
    df_file = info['file']
    df_pop = info['pop']
    
    info = df_name, path_to, df_key, df_file, df_pop
    
    return info


def _get_matrix_combos(geographic_units, non_geographies,
                       distance_type, representation):
    """generate all cost matrix combinations for calculation
    
    Parameters
    ----------
    geographic_units : list
        tracts, blocks groups, blocks?, populated blocks
    non_geographies : list
        incidents, fire stations
    distance_type : list
        distance / network allocation type
    representation : list
        centroid, pp2, va2n
    
    Returns
    -------
    combos : list
        all combinations of distance and observation type
    """
    
    combos = []
    
    #geog to same geo with same method
    for gu in geographic_units:
        for dt in distance_type:
            for rp in representation:
                combos.append([rp, dt, gu])
    
    #service to geo with same method
    for non_geo in non_geographies:
        for gu in geographic_units:
            for dt in distance_type:
                for rp in representation:
                    combos.append([rp, dt, non_geo, gu])
    
    return combos




