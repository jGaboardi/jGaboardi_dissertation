"""
"""

__author__ = 'James D. Gaboardi <jgaboardi@gmail.com>'

# standard library imports
import copy, os, time, warnings

# non-standard library imports
import geopandas as gpd
import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Point, MultiPoint
from shapely.geometry import LineString, MultiLineString
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import polygonize

# project library imports
from . import utils


def high_precision_id(c1, c2, c3):
    """mid level function -- vectorized function for generating a high
    precision object ID by combining the polygon ID, line segment ID,
    and dataframe index
    
    Parameters
    ----------
    c1 : cell in geopandas.GeoSeries
        GEOID
    c2 : cell in geopandas.GeoSeries
        line segment id
    c3 : cell in geopandas.GeoSeries
        dataframe index
    
    Returns
    -------
    idx : str
        combined GEOID_SegmID_index for pp2n or okabe
    """
    
    idx = '%s_%s_%s' % (c1, c2, c3)
    
    return idx


def pp2n(segms, polys, columns=None, offset_=None, remove_segm=None,
         restrict_col=None, drop_cols=None, xyid=None, geo_col=None,
         sid_name=None, poly_key=None, poly_pop=None, proj_init=None,
         line_pos=0.5, exhaustion_threshold=.9, pp2n_pointloc_increment=0.05,
         parallel_increment=0.5, pop_type=float, pp2n_id='pp2n_id',
         pp2n_col='pop_pp2n', mtfcc='MTFCC', len_segm='len_segm',
         len_tot='len_tot', pop_rat='pop_rat'):
    """top level function
    
       -- Population Polygon To Network (pp2n) --
          generate representative points near network segments
          which are assigned a proportional population weight
          based on segment length and total segment length of
          segments associated with census geographies.
    
    Parameters
    ----------
    segms : geopandas.GeoDataFrame
        road network segments
    polys : geopandas.GeoDataFrame
        census geoegraphies
    columns : list
        names of columns to create in new pp2n dataframe
    offset_ : float
        distance of pp2n parallel offset from the line segment.
        Default is None.
    line_pos : float
        normalized position along the line segment at which the initial
        attempt to create the pp2n point should be conducted.
        Default is 0.5.
    exhaustion_threshold : float
        normalized position along the line segment at which the
        algorithm should terminate if no feasible points are found.
        Default is 0.9.
    pp2n_pointloc_increment : float
        if a feasible pp2n point cannot be located at the `line_pos`
        try postive and negative increments of this value until either
        the line is exhausted or feasible locations are found.
        Default is 0.05.
    parallel_increment:
        decrease by `parallel_increment` when generating a parallel
        offset that results in a LinearRing or LINESTRING EMPTY.
        Default is 0.5.
    remove_segm : list
        segment MTFCC types to exclude from pp2n generation.
    restrict_col : str
        segment label to restrict (road segemnt mtfcc).
    drop_cols : list
        column names to drop following spatial join of the pp2n points
        and census geographies.
    pop_type : {float, int}
        data type of the generated pp2n population. Default is float.
    sid_name : str
        segment column name. Default is None.
    xyid : str
        combined x-coord + y-coords string ID. Default is None.
    proj_init : int
        intial coordinate reference system. default is None.
    geo_col : str
        geometry column name. Default is None.
    poly_key : str
        polygon dataframe key column name. Default is None.
    poly_pop : str
        polygon dataframe population column name. Default is None.
    pp2n_id : str
        pp2n ID column name. Default is 'pp2n_id'.
    pp2n_col : str
        pp2n population column name. Default is 'pp2n_col'.
    mtfcc : str
        MTFCC dataframe column name. Default is 'MTFCC'.
    len_segm : str
        segment length column name. Default is 'len_segm'.
    len_tot : str
        total associated length column name. Default is 'len_tot'.
    pop_rat : str
        population ratio column name. Default is 'pop_rat'.
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        generated, weighted pp2n points
    """
    
    # list of pp2n points
    pp2n_pts = create_pp2n(segms, polys,
                           offset_=offset_, line_pos=line_pos,
                           parallel_increment=parallel_increment,
                           exhaustion_threshold=exhaustion_threshold,
                           pp2n_pointloc_increment=pp2n_pointloc_increment,
                           remove=remove_segm, geo_col=geo_col,
                           mtfcc=mtfcc, sid_name=sid_name)
    
    # create pp2n dataframe
    df = create_pp2n_df(pp2n_pts, polys, segms, columns, crs=proj_init,
                        remove=remove_segm, drop=drop_cols, sid_name=sid_name,
                        pp2n_col=pp2n_col, len_segm=len_segm, len_tot=len_tot,
                        pop_rat=pop_rat, mtfcc=mtfcc, geo_col=geo_col,
                        restrict_col=restrict_col)
    
    # add unique id
    df[pp2n_id] = np.vectorize(high_precision_id)\
                              (df[poly_key], df[sid_name], df.index)
    
    # proportionalize population
    df = div_pop(df, pp2n_col=pp2n_col, poly_key=poly_key,
                 len_tot=len_tot, pop_rat=pop_rat, len_segm=len_segm,
                 poly_pop=poly_pop, as_type=pop_type)
    
    # due rounding precision during the previous step, some
    # pp2n point will be worth 0 population (i.e. population
    # ratio is 0.002, so population will be 0.0); drop these
    # cases as population total will be valid.
    df = df[df[pp2n_col] > 0.0]
    df.reset_index(inplace=True, drop=True)
    
    # set xyID
    node2xyid = utils.generate_xyid(df=df, geom_type='node', geo_col=geo_col)
    df = utils.fill_frame(df, idx='index', col=xyid, data=node2xyid)
    
    return df


def create_pp2n(segms, polys, offset_=None, line_pos=None,
                exhaustion_threshold=None, pp2n_pointloc_increment=None,
                parallel_increment=None, remove=None, geo_col=None,
                mtfcc=None, sid_name=None):
    """mid level function -- create a list of feasible pp2n points
    in the form -- [[idx1, point1], [idx1, point2], ...]
    
    Parameters
    ----------
    segms : see pp2n()
    polys : see pp2n()
    offset_ : see pp2n()
    line_pos : see pp2n()
    exhaustion_threshold : see pp2n()
    pp2n_pointloc_increment : see pp2n()
    parallel_increment : see pp2n()
    remove : see `remove_segm` in pp2n()
    geo_col : see pp2n()
    mtfcc : see pp2n()
    sid_name : see pp2n()
    
    Returns
    -------
    pp2ns : list
        generated pp2n points and associated segment ids
    """
    
    # segment count
    n_segms = segms.shape[0]
    
    # unary union of all census geographies in the study area
    union_polys = create_valid_unary(polys.unary_union)
    
    # list to store pp2n points
    pp2ns = []
    
    # iterate over each segment in the set
    for idx in range(n_segms):
        
        # skip if not parallel offset should be attempted
        if segms.loc[idx, mtfcc] in remove:
            continue
        
        # line of interest
        loi = segms.loc[idx, geo_col]
        
        # create a buffered offset from line-of-interest enveope
        loi_buffer = loi.envelope.buffer(offset_)
        
        # then extract the intersection
        local_poly = union_polys.intersection(loi_buffer)
        
        # prevent shapely from rearranging the line coords
        # (as in the case when a creating a new linear ring)
        is_ring = False
        
        if loi.is_ring:
            is_ring = True
        
        # create two parallel offsets of the `loi`
        offset_lines = gen_offsets(loi, offset_, for_ring=is_ring,
                                   parallel_increment=parallel_increment)
        
        if not offset_lines:
            seg_id = segms.loc[idx, sid_name]
            warn_msg = '\t\tno parallel lines generated for %s' % seg_id
            warnings.warn(warn_msg)
            continue
        
        determining_pp2n = True
        # iterate over the two offset lines until (1) two points
        # have been generated that intersect with the unioned
        # geographies study area; or (2) the line length is
        # exhausted.
        while determining_pp2n:
            segm_pp2n_pts = []
            
            for ofs_idx, ofl in enumerate(offset_lines):
                fixing_loc = True
                pp2n_pointloc = line_pos
                
                while fixing_loc:
                    is_midpoint = False
                    
                    # if the generated pp2n point is at the center
                    # of the line, record it and proceed
                    if pp2n_pointloc == line_pos:
                        
                        # interpolate pp2n point along the line
                        pp2n_point = ofl.interpolate(pp2n_pointloc,
                                                     normalized=True)
                        
                        # if the point is within the single, unioned
                        # polygons object record it
                        if pp2n_point.intersects(local_poly):
                            segm_pp2n_pts.append(pp2n_point)
                            is_midpoint = True
                            
                            # once two pp2n points have been located
                            # break out of the while loop
                            if len(segm_pp2n_pts) == 2:
                                determining_pp2n = False
                                
                            # once a viable point has been located
                            # break out of the location fixer
                            fixing_loc = False
                            
                    # if the generated pp2n point is not at the
                    # center of the line, iterate in +/- .05
                    # increments until an intersecting point is
                    # found or the line is exhausted
                    if not is_midpoint:
                        while pp2n_pointloc < exhaustion_threshold\
                        and fixing_loc:
                            
                            pp2n_pointloc += pp2n_pointloc_increment
                            increments = [pp2n_pointloc*-1.,
                                          pp2n_pointloc]
                            
                            for incr in increments:
                                
                                # interpolate pp2n point along the line
                                pp2n_point = ofl.interpolate(incr,
                                                             normalized=True)
                                
                                if pp2n_point.intersects(local_poly):
                                    segm_pp2n_pts.append(pp2n_point)
                                    
                                    # once a viable point has been
                                    # located break out of the
                                    # location fixer
                                    fixing_loc = False
                                    
                        # once two pp2n points have been located
                        # break out of the while loop
                        if len(segm_pp2n_pts) == 2:
                            determining_pp2n = False
                        
                        # if not viable point has been located
                        # after exhausting the line, break out
                        # of the location fixer
                        fixing_loc = False
                        
            # break out of loop when either (a) pp2n points have
            # been located; or (b) pp2n points have been deteremined
            # to not be feasible
            determining_pp2n = False
        
        # add pp2n points to list unless none were located
        if segm_pp2n_pts:
            pp2ns.extend([[idx, pp2n_pt] for pp2n_pt in segm_pp2n_pts])
    
    return pp2ns
    
    
def create_valid_unary(unary, nominal_buffer=0.):
    """if polygon is no valid make valid
    this works if the object is a Polygon or MultiPolygon
    
    Parameters
    ----------
    unary : {shapely.Polygon, shapely.MultiPolygon}
        may or may not consist of valid geometries
    nominal_buffer : {float, int}
        0. or small buffer to perform in order to revalid
        polygon geometries
    
    Returns
    -------
    unary : {shapely.Polygon, shapely.MultiPolygon}
        single valid geometries or multiple valid geometries
    """
    
    # if the geoemtry(ies) are already valid do nothing
    if not unary.is_valid:
        try:
            
            # if there are multiple polygons this will pass, else
            # as TypeError will be thrown because
            # 'Polygons have no length'
            len(unary)
            
            # convert MultiPolygon to list of polygons
            unary_polys = list(unary)
            valid_polys = []
            
            for poly in unary_polys:
                
                # if not valid geometry, make valid by
                # performing a buffer
                if not poly.is_valid:
                    poly = poly.buffer(nominal_buffer)
                
                # add the valid geometry to list
                if poly.is_valid:
                    valid_polys.append(poly)
                
                # if geometry still not valid we have a problem
                else:
                    raise ValueError('Polygon still not valid after buffer.')
            
            # cast back as MultiPolygon
            unary = MultiPolygon(valid_polys)
        
        except TypeError:
            unary = unary.buffer(nominal_buffer)
    
    return unary


def gen_offsets(loi_, offset_ , parallel_increment=None, for_ring=False):
    """low level function for generating parallel offset lines
    
    Parameters
    ----------
    loi_ : shapely.geometry.LineString
        network line segment in question
    offset_ : see pp2n(offset_)
    parallel_increment : see pp2n()
    for_ring : bool
        flag if segment is a LinearRing [True] or not [False].
        Default is False.
    
    Returns
    -------
    offset_lines : list
        two parallel line offset by `offset_` meters
    """
    
    # standard 'non-ring' line case
    if not for_ring:
        
        offset_lines = _gen_offsets(loi_,
                                    offset_,
                                    parallel_increment)
    
    # linear ring case
    else:
        
        fixing_ring_inset = for_ring
        while fixing_ring_inset:
            
            # generate parallel offset lines
            loi_ = LineString(loi_.coords[1:-1])
            offset_lines = _gen_offsets(loi_,
                                        offset_,
                                        parallel_increment)
            still_ring = False
            
            # check to see if the inner offset
            # resulted in a LinearRing
            for ofl in offset_lines:
                
                # if it did decrease the offset
                # by `parallel_increment`
                if ofl.is_ring:
                    offset_ -= parallel_increment
                    still_ring = True
                
            # segment no long a ring
            # break out of loop
            if not still_ring:
                fixing_ring_inset = False
    
    return offset_lines


def _gen_offsets(loi_, initial_offset, parallel_increment):
    """helper function for `gen_offsets()`
    
    Parameters
    ----------
    loi_ : see gen_offsets()
    initial_offset : see pp2n(offset_)
    parallel_increment : ee pp2n(parallel_increment)
    
    Returns
    -------
    offset_lines_ : list
        two parallel line offset by `initial_offset` meters or by
        `initial_offset` - (`parallel_increment` * INCR) when
        `initial_offset` results in an infeasible line; where INCR in
        the number of increments
    """
    
    offset_lines_ = []
    
    sides = ['right', 'left']
    
    for side in sides:
        complete = False
        offset = initial_offset
        
        while not complete:
            offset_loi = loi_.parallel_offset(offset, side=side)
            
            try:
                if len(offset_loi.coords[:]) > 1:
                    offset_lines_.append(offset_loi)
                    complete = True
                else:
                    offset -= parallel_increment
            
            except NotImplementedError:
                if type(offset_loi) == MultiLineString:
                    offset -= parallel_increment
                    if offset <= 0:
                        complete = True
                else:
                    raise TypeError(type(offset_loi), 'not supported')
    
    return offset_lines_


def create_pp2n_df(pts, polys, segms, cols, crs=None, remove=None, drop=None,
                   pp2n_col=None, geo_col=None, sid_name=None, len_segm=None,
                   len_tot=None, pop_rat=None, mtfcc=None, restrict_col=None):
    """mid level function -- instantiate a GeoDataFrame
    for the pp2n points
    
    Parameters
    ----------
    crs : dict
        coordinate reference system of the polygon dataframe
    pts : see pp2n()
    polys : see pp2n()
    segms : see pp2n()
    cols : see pp2n(columns)
    remove : see pp2n(remove_segm)
    drop : see pp2n(drop_cols)
    pp2n_col : see pp2n()
    restrict_col : see pp2n()
    sid_name : see pp2n()
    geo_col : see pp2n()
    len_segm : see pp2n()
    len_tot : see pp2n()
    pop_rat : see pp2n()
    mtfcc : see pp2n()
    
    Returns
    -------
    df : list
        generated pp2n points and associated segment ids
    """
    
    # columns for dataframe generation
    cols = cols + [pp2n_col, len_segm, len_tot, pop_rat]
    
    # generate dataframe of pp2n points
    empty_shell = np.empty((len(pts), len(cols)))
    df = gpd.GeoDataFrame(empty_shell, columns=cols, crs=crs,
                          geometry=[geom for (idx,geom) in pts])
    df = utils.set_crs(df, proj_init=crs)
    
    # set associated segment ids
    df[sid_name] = [idx for (idx,geom) in pts]
    
    # set associated segment lengths
    df[len_segm] = [segms.loc[idx, geo_col].length for idx in df[sid_name]]
    
    # record associated segment mtfcc
    df[restrict_col] = [segms.loc[idx, mtfcc] for idx in df[sid_name]]
    
    # perform a spatial join
    polys = utils.set_crs(polys, proj_init=crs)
    df = gpd.sjoin(df, polys)
    df = utils.set_crs(df, proj_init=crs)
    
    # drop irrelevant columns
    if drop:
        df.drop(drop, axis=1, inplace=True)
    
    return df


def div_pop(df, as_type=None, pp2n_col=None, poly_key=None, len_segm=None,
            poly_pop=None, len_tot=None, pop_rat=None):
    """ mid level function -- assign a ratio of the total population
    to each pp2n point within a census geography based on the ratio of
    the associated line segment length to the total length of line
    segments associated with the census geography
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        pp2n dataframe
    as_type : {int, float}
        data type of the generated pp2n population. Default is None.
    pp2n_col : see pp2n()
    poly_key : see pp2n()
    poly_pop : see pp2n()
    len_segm : see pp2n()
    len_tot : see pp2n()
    pop_rat : see pp2n()
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        updated pp2n dataframe
    """
    
    # for each group of pp2n points within
    # a single census geography
    for idx in df[poly_key].unique():
        
        # subset the dataframe
        pre_subset = df[df[poly_key] == idx].copy()
        
        # and calculate the total length of segments
        # associated with the census geography
        total_length = pre_subset[len_segm].sum()
        
        # then calculate the ratio of length and set the
        # calculate the population from the ratio
        for ss_idx in pre_subset.index:
            
            # set total length associated with geometry
            df.loc[ss_idx, len_tot] = total_length
            
            # set population ratio of associated geometry
            df.loc[ss_idx, pop_rat] = pre_subset.loc[ss_idx, len_segm] \
                                      / df.loc[ss_idx, len_tot]
            
            # set population associated with pp2n point
            df.loc[ss_idx, pp2n_col] = as_type(df.loc[ss_idx, pop_rat] \
                                                * df.loc[ss_idx, poly_pop])
        
        # check and adjust population as needed
        df = check_pop(df, idx, poly_key, poly_pop, pp2n_col)
    
    if as_type is int:
        df[pp2n_col] = df[pp2n_col].astype(as_type)
    
    # run final confirmation check for population
    _check_pop(df, idx, poly_key, poly_pop, pp2n_col)
    
    return df


def check_pop(df_, id_, poly_key, pop_col, pp2n_col,
              problems=None, round2=False):
    """low level function -- ensure the sum of the pp2n population
    is equal to the original population
    
    Parameters
    ----------
    df_ : geopandas.GeoDataFrame
        pp2n dataframe
    idx : int
        census geography ID
    poly_key : str
        census geography ID column name
    pop_col : str
        census geography population column name
    pp2n_col : see pp2n()
    problems : dict
        cases where the sum of ratio population does
        not equal the original total
    round2 : bool
        flag set for [True] if second round confirmation.
        Default is [False].
    
    Returns
    -------
    df_ : geopandas.GeoDataFrame
        adjusted (or not) pp2n dataframe
    problems : dict
        updated `problems` (only during second check)
    """
    
    subset = df_[df_[poly_key] == id_]
    pp2n_sum = subset[pp2n_col].sum()
    adjuster = subset.index[0]
    orig_sum = subset[pop_col][adjuster]
    
    if pp2n_sum != orig_sum:
        
        if not round2:
            diff = pp2n_sum - orig_sum
            actual = df_.loc[adjuster, pp2n_col]
            df_.loc[adjuster, pp2n_col] = actual - diff
        
        elif round2:
            problems[id_] = {'pp2n_sum':pp2n_sum,
                             'orig_sum':orig_sum}
            
    if not round2:
        return df_
    
    if round2:
        return problems


def _check_pop(df_, id_, poly_key, pop_col, pp2n_col,
               problems=None, round2=True):
    """ low level function -- ensure the sum of the pp2n population
    is equal to the original population
    
    Parameters
    ----------
    df_ : see check_pop()
    id_ : see check_pop()
    poly_key : see check_pop()
    pop_col : see check_pop()
    pp2n_col : see check_pop()
    problems : see check_pop()
    round2 : see check_pop()
    """
    
    problems = {}
    for idx in df_[poly_key].unique():
        problems = check_pop(df_, id_, poly_key, pop_col, pp2n_col,
                             problems=problems, round2=round2)
    
    if problems:
        warnings.warn('\n\npp2n sum is not equal to original population sum.')
        warnings.warn(problems.__str__())


def va2n(area, alloc_dir, geogs, net_segms, net_nodes, remove_segm=None,
         mtfcc=None, cols_in=None, geo_col=None, sid_name=None, inter=None,
         xyid=None, voronoi_diagnostic=True, clip_by=None, bounds=None,
         poly_key=None, poly_pop=None, va2n_polys=None, proj_init=None,
         initial_rho=10, incremental_rho=10, offset_param=2., area_thresh=.8,
         symdiff_buff=0.00001, halt_at=1000, va2n_id='va2n_id',
         va2n_col='pop_va2n', geo_area='geo_area', rat_area='rat_area',
         file_type='.shp'):
    """ Top level function for vector area-to-network data conversion method
    of allocating polygon weights to a network.
    
    CITE MOIOKA?
    
    Parameters
    ----------
    area : str
        study area within county. Default is None.
    alloc_dir : str
        path to `allocation` data. Default is None.
    geogs : geopandas.GeoDataFrame
        census geographies
    net_segms : geopandas.GeoDataFrame
        road network segments
    net_nodes : geopandas.GeoDataFrame
        road network nodes (articulation points)
    remove_segm : list
        no population on roads that should never have population.
    mtfcc : str
        MAF/TIGER Feature Class code. Default is None.
    cols_in : list
        columns to keep after spatial join / dissolve. Default is None.
    geo_col : str
        geometry column name. Default is None.
    sid_name : str
        segment column name. Default is None.
    inter : str
        file path to intermediary data. Default is None.
    xyid : str
        combined x-coord + y-coords string ID. Default is None.
    voronoi_diagnostic: bool
        print iteration diagnostics. Default is True.
    clip_by : shapely.geometry.[Multi]Polygon
        final boundary for trimming voronoi cells. Default is None.
    bounds : shapely.geometry.Polygon
        outer boundary to trim reaching Voronoi cells. Default is None.
    poly_key : str
        polygon dataframe key column name. Default is None.
    poly_pop : str
        polygon dataframe population column name. Default is None.
    va2n_polys : str
        derivation from census geography file name. Default is None.
    proj_init : int
        intial coordinate reference system. default is None.
    initial_rho : int
        intial count of points to generate along line segments.
        Default is 10.
    incremental_rho : int
        increase `initial_rho` by this number for each failed iteration
        of LVD generation, where a failed iteration is defined as an
        unequal number of line segment and LVD cells. Default is 10.
    area_thresh : {int, float}
        ratio of the calculated symmetric differnce areas of the LVD
        polygon cells to the area of the segments envelope. The
        `area_ratio` must be greater than `area_thresh` for the
        algorithm to terminate successfully. Default is .95.
    symdiff_buff : {int, float}
        value to buffer symmetrical difference overlay. Default is .5.
    halt_at : int
        break iteration. Default is 1000 (rho).
    va2n_id : str
        va2n ID column name. Default is 'va2n_id'.
    va2n_col : str
        va2n population column name. Default is 'va2n_col'.
        portion of population.
    geo_area : str
        label for column of area from census geographies subset.
        Default is 'geo_area'.
    rat_area : str 
        label for column of area ratio from census geographies subset.
        Default is 'rat_area'.
    file_type : str
        file extension. Default is '.shp'.
    
    Returns
    -------
    va2n_points_df : geopandas.GeoDataFrame
        one point at the mid point of each line segment that has
        a population associated with it.
    """
    
    # -- Line Voronoi Diagram
    lvd_file = '%slvd_%s_%s%s' % (alloc_dir, initial_rho, area, file_type)
    # if the file already exists skip
    lvd_file_exists = os.path.exists(lvd_file)
    if not lvd_file_exists:
    
        # create line voronoi diagram
        lvd = faux_lvd(area, net_segms, net_nodes, remove_segm=remove_segm,
                       mtfcc=mtfcc, bounds=bounds, clip_by=clip_by,
                       initial_rho=initial_rho, offset_param=offset_param,
                       incremental_rho=incremental_rho, id_col=sid_name,
                       cols_in=cols_in, symdiff_buff=symdiff_buff, 
                       diagnostic=voronoi_diagnostic, area_thresh=area_thresh,
                       geo_col=geo_col, inter=inter, proj_init=proj_init,
                       file_type=file_type)
        lvd.to_file(lvd_file)
    else:
        lvd = gpd.read_file(lvd_file)
        lvd = utils.set_crs(lvd, proj_init=proj_init)
    
    
    # -- va2n Polygons
    va2n_poly_file = '%s%s%s' % (inter, va2n_polys, file_type)
    # if the file already exists skip
    va2n_poly_file_exists = os.path.exists(va2n_poly_file)
    if not va2n_poly_file_exists:
        
        # overlay -- Union of census geographies and LVD cells 
        va2n_polys_df = gpd.overlay(geogs, lvd, how='union')
        va2n_polys_df = utils.set_crs(va2n_polys_df, proj_init=proj_init)
        # drop NaNs
        va2n_polys_df.dropna(inplace=True)
        # Convert IDs back to integer
        va2n_polys_df[poly_key] = va2n_polys_df[poly_key].astype(int)
        va2n_polys_df[sid_name] = va2n_polys_df[sid_name].astype(int)
        
        # create new empty column to house new va2n area
        va2n_polys_df[geo_area] = np.zeros
        
        # fetch area for each original census geography
        for unique_geoid in va2n_polys_df[poly_key].unique():
            subset = va2n_polys_df[va2n_polys_df[poly_key] == unique_geoid]
            subset_area = subset.area.sum()
            va2n_polys_df.loc[subset.index, geo_area] = subset_area
        
        # calculate va2n_polys / original_poly area ratio
        va2n_polys_df[rat_area] = va2n_polys_df.area\
                                  / va2n_polys_df[geo_area]
        # set va2n population proportion based on area ratio
        va2n_polys_df[va2n_col] = va2n_polys_df[rat_area]\
                                  * va2n_polys_df[poly_pop]
        
        # generate high-precision ID (polygon, segment, index)
        va2n_polys_df[va2n_id] = np.vectorize(high_precision_id)\
                                             (va2n_polys_df[poly_key],
                                              va2n_polys_df[sid_name],
                                              va2n_polys_df.index)
        va2n_polys_df.to_file(va2n_poly_file)
    
    else:
        va2n_polys_df = gpd.read_file(va2n_poly_file)
        va2n_polys_df = utils.set_crs(va2n_polys_df, proj_init=proj_init)
    
    del lvd
    
    
    # -- va2n Points
    va2n_point_file = '%sva2n_%s_%s%s' % (alloc_dir, initial_rho,
                                           area, file_type)
    # if the file already exists skip
    va2n_point_file_exists = os.path.exists(va2n_point_file)
    if not va2n_point_file_exists:
        
        # va2n point dataframe
        va2n_points_df = gpd.GeoDataFrame(crs=net_segms.crs)
        
        # segments unvisited set
        va2n_poly_seg_ids = set(va2n_polys_df[sid_name])
        
        while va2n_poly_seg_ids:
            
            # iterate over network segments
            for sid in net_segms[sid_name]:
                
                # evaluate only if the segment is associated
                # with an `va2n polygon`
                if sid in va2n_poly_seg_ids:
                    
                    # extract population of polygons
                    # associated with the segment
                    pop_bool = (va2n_polys_df[sid_name] == sid)
                    _pop_va2n_ = va2n_polys_df.loc[pop_bool, va2n_col]
                    _pop_va2n_ = _pop_va2n_.squeeze()
                    
                    # convert datatype and sum if more than one value
                    if type(_pop_va2n_) == float:
                        pass
                    elif type(_pop_va2n_) == str:
                        _pop_va2n_ = float(_pop_va2n_)
                    else:
                        _pop_va2n_ = _pop_va2n_.astype(float).sum()
                    
                    # extract IDs of polygons
                    # associated with the segment
                    ids_bool = (va2n_polys_df[sid_name] == sid)
                    _va2n_ids = list(va2n_polys_df.loc[ids_bool, va2n_id])
                    
                    # create new unique id based on the `okabe ids`
                    if len(_va2n_ids) > 1:
                        _va2n_ids = '-'.join(_va2n_ids)
                    else:
                        _va2n_ids = _va2n_ids[0]
                    
                    # dataframe index of segment
                    seg_df_idx = (net_segms[sid_name] == sid)
                    # segment geometry
                    segment = net_segms.loc[seg_df_idx, geo_col].squeeze()
                    # segment midpoint
                    va2n_point = segment.interpolate(.5, normalized=True)
                    
                    # set up the record and column structure
                    record = [[sid, va2n_point, _va2n_ids, _pop_va2n_]]
                    columns = [sid_name, geo_col, va2n_id, va2n_col]
                    
                    # prepare a temporary dataframe for the points
                    temp_gdf = gpd.GeoDataFrame(record, columns=columns)
                    
                    # append temporary dataframe to full dataframe
                    va2n_points_df = va2n_points_df.append(temp_gdf)
                    
                    # mark the segments as visited
                    va2n_poly_seg_ids.remove(sid)
        
        va2n_points_df.reset_index(inplace=True, drop=True)
        va2n_points_df = utils.set_crs(va2n_points_df, proj_init=proj_init)
        # set xyID
        node2xyid = utils.generate_xyid(df=va2n_points_df, geom_type='node',
                                        geo_col=geo_col)
        va2n_points_df = utils.fill_frame(va2n_points_df, idx='index',
                                          col=xyid, data=node2xyid)
    
    return va2n_points_df


def faux_lvd(area, segms, nodes, remove_segm=None, mtfcc=None, inter=None,
             geo_col=None, offset_param=None, bounds=None, clip_by=None,
             id_col=None, cols_in=None, diagnostic=True, initial_rho=None,
             incremental_rho=None, halt_at=None, symdiff_buff=None,
             area_thresh=None, proj_init=None, file_type=None):
    """Generate a faux Voronoi diagram from line segments (LVD).
    Voronoi cells are generated by creating incremental points along
    each line segment, thereby representing the line segment as a dense
    chain of points. With the network as this dense chain of points,
    generate a Voronoi diagram from points(PVD). Next, perform a spatial
    join on with the PVD and the network line segments. Finally,
    dissolve the PVD by the associated network segment ID attribute,
    resulting in LVD polygon cells. Following LVD creation, compare the
    number of LVD polygon cells to the initial segment count. If these
    numbers are equal return the LVD dataframe. If the numbers are not
    equal (a failed iteration), increase the density (rho) incrementally
    and perform another itertion of LVD creation.
    
    Parameters
    ----------
    area : str
        study area within county. Default is None.
    segms : geopandas.GeoDataFrame
        network line segments
    nodes : geopandas.GeoDataFrame
        network articulation nodes
    remove_segm : list
        no population on roads that should never have population
    mtfcc : str
        MAF/TIGER Feature Class code. Default is None.
    inter : str
        file path to intermediary data. Default is None.
    geo_col : str
        geometry column name. Default is None.
    offset_param : {int, float}
        distance from articulation point to begin and end
        segment densification. Default is 0.5 (meters).
    bounds : shapely.geometry.Polygon
        outer boundary to trim reaching Voronoi cells. Default is None.
    clip_by : shapely.geometry.[Multi]Polygon
        final boundary for trimming voronoi cells. Default is None.
    id_col : str
        segment dataframe ID column name. Default is None.
    cols_in : list
        columns to keep after spatial join / dissolve. Default is None.
    diagnostic : bool
        print iteration diagnostics. Default is True.
    initial_rho : int
        intial count of points to generate along line segments.
        Default is 10.
    incremental_rho : int
        increase `initial_rho` by this number for each failed iteration
        of LVD generation, where a failed iteration is defined as an
        unequal number of line segment and LVD cells. Default is 10.
    halt_at : int
        break iteration. Default is 1000 (rho).
    symdiff_buff : {int, float}
        value to buffer symmetrical difference overlay. Default is .5.
    area_thresh : {int, float}
        ratio of the calculated symmetric differnce areas of the LVD
        polygon cells to the area of the segments envelope. The
        `area_ratio` must be greater than `area_thresh` for the
        algorithm to terminate successfully. Default is .95.
    proj_init : int
        intial coordinate reference system. default is None.
    file_type : str
        file extension. Default is None.
    
    Raises
    ------
    RuntimeError
        Raised when a feasible solution is not found after a specified
        number of iterations/when the `rho` limit is reached.
    
    Returns
    -------
    vor_df : geopandas.GeoDataFrame
        faux Voronoi polygons generated from line segments
    """
    
    # create a outer boundary to clip the portions of
    # Voronoi polygon cells that either extend into infinity or are
    # out of the study area.
    if bounds is None:
        bounds =  MultiLineString(list(segms[geo_col])).envelope
    
    # inital total area of segments envelope
    if isinstance(clip_by, gpd.GeoDataFrame):
        init_area = clip_by.area.squeeze()
    else:
        init_area = bounds.area.squeeze()
    
    # count of line segments in the network
    if remove_segm and mtfcc:
        segment_count = segms[~segms[mtfcc].isin(remove_segm)].shape[0]
    else:
        segment_count = segms.shape[0]
    
    equal_count = False
    algo_complete = False
    iteration = 0
    
    while not algo_complete:
        
        # set rho
        rho = initial_rho + incremental_rho * iteration
        
        
        # Densified Line Segment Vertices
        vtxs_file = '%svtxs_df_%s_%s%s' % (inter, area, iteration, file_type)
        # if the file already exists skip
        vtxs_file_exists = os.path.exists(vtxs_file)
        if not vtxs_file_exists:
            
            # generate a dense point chain
            # representation of the line segments
            densified_vtxs = dense_vertices(segms,
                                            offset_param,
                                            rho=rho,
                                            id_col=id_col)
            
            # instantiate a dataframe for the dense points
            vtxs_df = gpd.GeoDataFrame(densified_vtxs,
                                       columns=[id_col, geo_col])
            vtxs_df = utils.set_crs(vtxs_df, proj_init=proj_init)
            if inter:
                vtxs_df.to_file(vtxs_file)
        else:
            vtxs_df = gpd.read_file(vtxs_file)
            vtxs_df = utils.set_crs(vtxs_df, proj_init=proj_init)
        
        
        # Point Voronoi Diagram Polygons
        pvd_polys_file = '%spvd_polys_df_%s_%s%s' % (inter, area,
                                                     iteration, file_type)
        # if the file already exists skip
        pvd_polys_file_exists = os.path.exists(pvd_polys_file)
        if not pvd_polys_file_exists:
            
            # generate point Voronoi diagram from dense point chains
            pvd = Voronoi([[p.x, p.y] for p in vtxs_df[geo_col]])
            
            # create dataframe of polygonized voronoi cell ridge edges
            pvd_polys_df = polygonize_ridges(pvd, bounds)
            pvd_polys_df = utils.set_crs(pvd_polys_df, proj_init=proj_init)
            del pvd
            
            if inter:
                pvd_polys_df.to_file(pvd_polys_file)
        else:
            pvd_polys_df = gpd.read_file(pvd_polys_file)
            pvd_polys_df = utils.set_crs(pvd_polys_df, proj_init=proj_init)
        
        del vtxs_df
        
        
        # Line Voronoi Diagram Polygons -- All
        lvd_polys_file = '%sall_lvd_polys_df_%s_%s%s' % (inter, area,
                                                         iteration, file_type)
        # if the file already exists skip
        lvd_polys_file_exists = os.path.exists(lvd_polys_file)
        if not lvd_polys_file_exists:
            
            # generate buffer of nodes to match the `offset_param`
            # this is performed to handle precision issues
            nodes_buffer = nodes.buffer(offset_param)
            nodes_buffer = gpd.GeoDataFrame(geometry=nodes_buffer)
            nodes_buffer = utils.set_crs(nodes_buffer, proj_init=proj_init)
            
            # spatial join and dissolve
            lvd_polys_df = spjoin_dissolve(pvd_polys_df, segms,
                                           nodes_buffer, id_col, cols_in)
            lvd_polys_df = utils.set_crs(lvd_polys_df, proj_init=proj_init)
            
            if inter:
                lvd_polys_df.to_file(lvd_polys_file)
        else:
            lvd_polys_df = gpd.read_file(lvd_polys_file)
            lvd_polys_df = utils.set_crs(lvd_polys_df, proj_init=proj_init)
        
        del pvd_polys_df
        
        # set break condition for `Test_Grid_Leon_FL`
        if area == 'Test_Grid_Leon_FL':
            algo_complete = True
            return lvd_polys_df
        
        
        # Line Voronoi Diagram Polygons -- Clipped and Cleaned
        lvd_polys_file = '%scac_lvd_polys_df_%s_%s%s' % (inter, area,
                                                         iteration, file_type)
        # if the file already exists skip
        lvd_polys_file_exists = os.path.exists(lvd_polys_file)
        if not lvd_polys_file_exists:
            
            # intersections
            inter_gdf = gpd.GeoDataFrame()
            inter_gdf = utils.set_crs(inter_gdf, proj_init=proj_init)
            # find intersecting and overlapping polygons
            lvd_polys_intersections = get_intersections(lvd_polys_df,
                                                        bounds,
                                                        inter_gdf)
            lvd_polys_intersections = utils.set_crs(lvd_polys_intersections,
                                                    proj_init=proj_init)
            
            # calculate the symmetric difference of
            # `lvd_polys_intersections` and then keep only the objects
            # produced from the original geoms
            # calculate buffer around the intersections objects
            df_buff = lvd_polys_intersections.buffer(symdiff_buff)
            df_buff = gpd.GeoDataFrame(geometry=df_buff)
            df_buff = utils.set_crs(df_buff, proj_init=proj_init)
            lvd_polys_df = get_symdiff(lvd_polys_df, df_buff, id_col)
            lvd_polys_df = utils.set_crs(lvd_polys_df, proj_init=proj_init)
            
            # remove all unnecessary columns
            clip_by = clip_by[clip_by.columns[clip_by.columns == geo_col]]
            # perform clip
            lvd_polys_df = gpd.overlay(lvd_polys_df, clip_by,
                                       how='intersection')
            lvd_polys_df = utils.set_crs(lvd_polys_df, proj_init=proj_init)
            
            if inter:
                lvd_polys_df.to_file(lvd_polys_file)
        else:
            lvd_polys_df = gpd.read_file(lvd_polys_file)
            lvd_polys_df = utils.set_crs(lvd_polys_df, proj_init=proj_init)
        
        
        # total area of symmetric difference LVD poltgon cells
        post_area = lvd_polys_df.area.sum()
        # ratio of the inital evelope area to the LVD area
        area_ratio = post_area / init_area
        
        
        # Line Voronoi Diagram Polygons -- Only Permissible MTFCC
        # remove polys than don't support population
        lvd_polys_df = lvd_polys_df[~lvd_polys_df[mtfcc].isin(remove_segm)]
        lvd_polys_file = '%ssub_lvd_polys_df_%s_%s%s' % (inter, area,
                                                         iteration, file_type)
        if inter:
            lvd_polys_df.to_file(lvd_polys_file)
        
        # iteration count and bools definition
        iteration += 1
        # number of polygons to compare with line segments
        lvd_polygon_count = lvd_polys_df.shape[0]
        equal_count_bool = segment_count == lvd_polygon_count
        
        # set break condition for `Test_Tract_Leon_FL`
        if area == 'Test_Tract_Leon_FL':
            equal_count_bool = True
            
        area_ratio_bool = area_ratio > area_thresh
        
        
        # Identify non-polygonal (including multi) geometries
        non_poly_ids = []
        for idx in lvd_polys_df.index:
            geom = lvd_polys_df.loc[idx, lvd_polys_df[geo_col].name]
            if type(geom) not in [Polygon, MultiPolygon]:
                only_polys = False
                non_poly_ids.append(lvd_polys_df.loc[idx, sid_name])
        
        # non-polygonal geometries found
        if diagnostic and non_poly_ids:
            msg1 = '\tTrouble segment IDs '
            msg2 = '(geometry other than Polygon/MultiPolygon): '
            print('%s%s%s' % (msg1, msg2, non_poly_ids))
        
        # algorithm successful -- break condition
        if equal_count_bool and area_ratio_bool:
            algo_complete = True
        
        # solution not found -- raise error
        if rho == halt_at:
            err_msg = 'Solution not found after %s iterations ' % iteration \
                      + 'with a final rho of %s.' % rho
            raise RuntimeError(err_msg)
        
        # solution not found -- start next iteration
        if diagnostic and not algo_complete:
            n_sp = (segment_count, lvd_polygon_count)
            r_incr = (rho, incremental_rho)
            areas = (rho, incremental_rho)
            print('After iteration %s: ' % iteration \
                  + 'segment count (%s) -- polygon count (%s) -- ' % n_sp \
                  + 'area ratio (%s)\n' % area_ratio \
                  + 'Increasing rho (%s) by %s.' % r_incr)
    
    return lvd_polys_df
    
    
def dense_vertices(segms, offset_param, rho=None, id_col=None):
    """Generate a dense chain of points for each network line segment.
    
    Parameters
    ----------
    segms : see `faux_voronoi_from_lines()`
    offset_param : see `faux_voronoi_from_lines()`
    rho : see `faux_voronoi_from_lines()`
    id_col : see `faux_voronoi_from_lines()`
    
    Returns
    -------
    dense_vtxs : list
        dense point chain representation of line segments tagged
        with the parent line segment ID in the form
        [[idx1, geom1], [idx1, geom2].
    """
    
    dense_vtxs = []
    
    # iterate over segments and create points along each segment at each
    # location in `increments` and tag with the ID of the line segment
    # on which they are located.
    for df_idx in segms.index:
        
        # fetch ID and geometry
        idx = segms.loc[df_idx, id_col]
        geom = segms.loc[df_idx, segms.geometry.name]
        
        # total length of segment
        tot_len = geom.length
        #
        increments = np.linspace(offset_param, tot_len-offset_param, rho)
        
        # create points based on increments along the segment
        for point in increments:
            dense_vtxs.append([idx, geom.interpolate(point)])
        
    return dense_vtxs


def polygonize_ridges(vd, bounds):
    """First create polygons from the voronoi cell ridge edges.
    Then trim the resulting polygons that are outside the study area.
    Finally return a GeoDataFrame of the voronoi diagram polygon cells.
    
    Parameters
    ----------
    vd : scipy.spatial.Voronoi
        voronoi diagram object
    bounds : shapely.geometry.Polygon
        outer boundary of network segments
    
    Returns
    -------
    vd_polys_df : geopandas.GeoDataFrame
        voronoi diagram cells
    """
    
    # create lines from PVD ridges 
    ridge_lines = [LineString(vd.vertices[line])
                   for line in vd.ridge_vertices
                   if -1 not in line]
    
    # polygonize ridge lines
    vd_polys = polygonize(ridge_lines)
    vd_polys_df = gpd.GeoDataFrame(geometry=list(vd_polys))
    
    return vd_polys_df

def spjoin_dissolve(pvd_polys_df, segms, nodes_buffer, id_col, cols_in):
    """(1) perform a spatial join on voronoi diagram cells and the
    network segments; and (2) and dissolve the cells on the segment id.
    
    Parameters
    ----------
    pvd_polys_df : geopandas.GeoDataFrame
    segms : see `faux_voronoi_from_lines()`
    nodes_buffer : geopandas.GeoDataFrame
        buffer of nodes to match the `offset_param`
    offset_param : see `faux_voronoi_from_lines()`
    id_col : see `faux_voronoi_from_lines()`
    cols_in : see `faux_voronoi_from_lines()`
    
    Returns
    -------
    lvd : geopandas.GeoDataFrame
        cells of a voronoi diagram generated from lines
    """
    
    # overlay of point-voronoi and node-buffer
    pvd_polys_df = gpd.overlay(pvd_polys_df,
                               nodes_buffer,
                               how='difference')
    
    # spatial join of point-voronoi and line semgents
    joined_df = gpd.sjoin(pvd_polys_df, segms)
    
    # drop some columns
    cols_in = [joined_df.geometry.name] + cols_in
    joined_df = joined_df[cols_in]
    
    # dissolve point-voronoi based on line segment id
    lvd_polys_df = joined_df.dissolve(id_col)
    lvd_polys_df[id_col] = lvd_polys_df.index
    
    return lvd_polys_df


def get_intersections(df, bounds, inter_gdf):
    """ iterate over all combinations of geometries produced by
    `spjoin_dissolve()` and isolate intersections.
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        result from spjoin_dissolve()
    bounds : geopandas.GeoDataFrame
        full bounding box of study area
    inter_gdf : geopandas.GeoDataFrame
        empty geopandas.GeoDataFrame to store intersections
    
    Returns
    -------
    inter_gdf : list
        intersection of intersecting geometries from `df`.
    """
    
    # Polygon, MultiPolygon intersections are not acceptable
    acceptable_intersections = [Point, MultiPoint,
                                LineString, MultiLineString]
    
    # create quadrat from bounds
    quad_polys = create_quadrat(bounds)
    
    # spatial index
    sindex = df.sindex
    
    # tracker for only evaluating pairs once
    visited_pairs = set()
    
    for poly in quad_polys:
        # buffer by the <1 micron dist to account for any space lost
        # in the quadrat cutting  otherwise may miss point(s) that
        # lay directly on quadrat line
        poly = poly.buffer(np.finfo(float).eps).buffer(0)
        
        # find approximate matches with r-tree, then precise
        # matches from those approximate ones
        possible_matches_index = list(sindex.intersection(poly.bounds))
        possible_matches = df.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(poly)]
        
        # intersecting/overlapping geometries
        overlaps = []
        for i, idx in enumerate(precise_matches.index):
            for j, jdx in enumerate(precise_matches.index):
                
                # record i,j pair
                pair = tuple(sorted((idx, jdx)))
                if idx == jdx:
                    continue
                if pair in visited_pairs:
                    continue
                
                # extract geometries
                i_geom = precise_matches.loc[idx, df.geometry.name]
                j_geom = precise_matches.loc[jdx, df.geometry.name]
                
                # determine the intersecting relationship
                if i_geom.intersects(j_geom):
                    relation = i_geom.intersection(j_geom)
                    if type(relation) not in acceptable_intersections:
                        overlaps.append(relation)
                
                # mark pair as evaluated
                visited_pairs.add(pair)
        
        # append the geometries to the initial dataframe
        inter_gdf = inter_gdf.append(gpd.GeoDataFrame(geometry=overlaps))
    
    return inter_gdf


def create_quadrat(bounds, grid_dim=10):
    """cut up a bounding polygon into quadrat
    
    Parameters
    ----------
    bounds : geopandas.GeoDataFrame
        full bounding box of study area
    grid_dim : int
        rows and columns count
    
    Returns
    -------
    polys : list
        shapely Polygon/Multipolygon objects as quadrats
        of bounding area
    """
    
    # Point, LineString, etc. intersections are not acceptable
    acceptable_intersections = [Polygon, MultiPolygon]
    
    # quad rows
    rows = np.linspace(bounds.bounds.miny,
                       bounds.bounds.maxy,
                       grid_dim)
    
    # quad columns
    cols = np.linspace(bounds.bounds.minx,
                       bounds.bounds.maxx,
                       grid_dim)
    
    # iterate over each quad to extract the intersection
    polys = []
    for ridx, row in enumerate(rows[:-1]):
        for cidx, col in enumerate(cols[:-1]):
            
            # lower left corner
            origin = (col, row)
            # lower right corner
            lr = (cols[cidx+1], row)
            # top right corner
            tr = (cols[cidx+1], rows[ridx+1])
            # top left corner
            tl = (col, rows[ridx+1])
            
            # polygonize quadrat
            poly = Polygon(LineString((origin, lr, tr, tl, origin)))
            
            # get intersection of quadrat and bounds area
            poly = poly.intersection(bounds.geometry.squeeze())
            
            # only keep polygons and multipolygons
            if type(poly) in acceptable_intersections:
                polys.append(poly)
    
    return polys


def get_symdiff(lvd, df_buff, id_col):
    """(1) perform a symmetric difference of intersecting geometries
    buffer and LVD; and (2) keep only geometries produced by
    original objects.
    
    Parameters
    ----------
    lvd : geopandas.GeoDataFrame
    df_buff : geopandas.GeoDataFrame
        bufferec product of `get_intersections()`
    id_col : see `faux_voronoi_from_lines()`
    
    Returns
    -------
    lvd : geopandas.GeoDataFrame
        non-intersecting LVD
    """
    
    # symmetric difference overlay to find overlaps
    lvd = gpd.overlay(lvd, df_buff, how='symmetric_difference')
    
    # drop objects not produced by original geometries
    lvd.dropna(inplace=True, axis=0)
    
    # recast the ids columns as integer
    lvd[id_col] = lvd[id_col].astype(int)
    
    return lvd
