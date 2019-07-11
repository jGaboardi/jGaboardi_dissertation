"""General utilities for the dissertation at1866 workflow
"""

# standard library imports
import copy, os, re, subprocess, sys, time

# non-standard library imports
import geopandas as gpd
import numpy as np
import pandas as pd
#import pysal as ps
import shapely


class ObjectSize:
    """retreive size of an object. Default set in GB.
    """
    
    
    def __init__(self, obj, single_object=False, units='GB', total_rounder=5):
        """retreive size of an object. Default in GB.
        1073741824.0005517 bytes in one gigabyte
        
        Parameters
        ----------
        obj : any python object
        single_object : bool
            get size of all attributes within object [True], or a single
            object/attribute [False]. Default is False.
        units : str
            unit to measure object size. Default is 'GB'.
        total_rounder : int
            decimal places to round the total object size. Default is 5.
        """
        
        self.time_snapshot = time.strftime('%Y-%m-%d %H:%M')
        self.obj = obj
        self.units = units
        self.total_rounder = total_rounder
        self.col_str = 'size_'+self.units
        if units == 'GB':
            self.unit_conversion = 1073741824.0005517
        
        
        def _attr_by_type(self, attrs_type=None):
            """get private and public attributes within an objects
            
            Parameters
            ----------
            attrs_type : str
                type of attribute requestes. Options are 'public' and
                'private'. Default is None.
            
            Returns
            -------
            attrs : list
                list of attribute in object
            """
            
            if attrs_type == 'public':
                attrs = [a for a in dir(self.obj) if not a.startswith('_')]
            
            if attrs_type == 'private':
                attrs = [a for a in dir(self.obj) if a.startswith('_')]
            
            else:
                attrs = [a for a in dir(self.obj)]
            
            return attrs, len(attrs)
        
        
        def _round_tot(self, tot):
            """round the toal object size to 'tot' decimal places
            
            Parameters
            ----------
            tot : int
                number of decimal places to round the total object size.
            
            Returns
            -------
            rounded float
            """
            
            return round(tot, self.total_rounder)
        
        
        def get_indiv_object_size(self, attr):
            """
            calculate the size of an individual object/attribute
            
            Parameters
            ----------
            attr : str
            
            Returns
            -------
            size (in self.unit_conversion) of the individual attribute
            """
            
            return sys.getsizeof(attr)/self.unit_conversion
        
        
        def get_size_frame(self, a):
            """create pandas dataframe with object information
            
            Parameters
            ----------
            a : list
                list of attributes in object
            
            Returns
            -------
            frame : pandas.DataFrame
                dataframe containing information about object attributes
            """
            
            # columns for dataframe
            cols = ['attr', 'attr_type', 'type', self.col_str]
            
            # instantiate dataframe
            frame = pd.DataFrame(np.empty((len(a), len(cols))), columns=cols)
            
            # fill dataframe with information for each attribute
            for idx, attr_str in enumerate(a):
                attr = getattr(self.obj, attr_str)
                if attr_str.startswith('_'):
                    attr_type = 'private'
                else:
                    attr_type = 'public'
                
                # calculate individual size
                size = get_indiv_object_size(self, attr)
                
                # add record -- name, type, string of type, size
                frame.iloc[idx] = attr_str, attr_type, str(type(attr)), size
            
            # set attribute name as the index
            frame.index = frame['attr']
            
            return frame
        
        if single_object:
            self.obj_size = get_indiv_object_size(self, self.obj)
        
        else:
            
            # get public attributes and public attribute count
            self.attrs_pub, n_attrs_pub = _attr_by_type(self,
                                                        attrs_type='public')
            
            # get private attributes and private attribute count
            self.attrs_prv, n_attrs_prv = _attr_by_type(self,
                                                        attrs_type='private')
            # get all attributes and all attribute count
            
            self.attrs_all, n_attrs_all = _attr_by_type(self)
        
            # size dataframes
            self.size_frame_pub = get_size_frame(self, self.attrs_pub)
            self.size_frame_prv = get_size_frame(self, self.attrs_prv)
            self.size_frame_all = get_size_frame(self, self.attrs_all)
            
            # size totals
            self.size_tot_pub = self.size_frame_pub[self.col_str].sum()
            self.size_tot_prv = self.size_frame_prv[self.col_str].sum()
            self.size_tot_all = self.size_frame_all[self.col_str].sum()
            
            # round total size
            if total_rounder:
                self.size_tot_pub = _round_tot(self, self.size_tot_pub)
                self.size_tot_prv = _round_tot(self, self.size_tot_prv)
                self.size_tot_all = _round_tot(self, self.size_tot_all)


def get_system_memory():
    """retreive and print system memory at point in time (gigabytes)
    """
    
    # Get process info
    ps = subprocess.Popen(['ps', '-caxm', '-orss,comm'],
                          stdout=subprocess.PIPE).communicate()[0].decode()
    vm = subprocess.Popen(['vm_stat'],
                          stdout=subprocess.PIPE).communicate()[0].decode()
    
    # Iterate processes
    processLines = ps.split('\n')
    sep = re.compile('[\s]+')
    rssTotal = 0 # kB
    
    for row in range(1,len(processLines)):
        rowText = processLines[row].strip()
        rowElements = sep.split(rowText)
        
        try:
            rss = float(rowElements[0]) * 1024
        except:
            
            rss = 0
        rssTotal += rss
    
    # Process vm_stat
    vmLines = vm.split('\n')
    sep = re.compile(':[\s]+')
    vmStats = {}
    
    for row in range(1,len(vmLines)-2):
        rowText = vmLines[row].strip()
        rowElements = sep.split(rowText)
        vmStats[(rowElements[0])] = int(rowElements[1].strip('\.')) * 4096
    
    print('\t\t=====================================')
    print('\t\tSystem Memory:')
    print('\t\t-------------')
    print('\t\tWired:\t\t\t%d GB'%(vmStats['Pages wired down']/1024/1024/1024))
    print('\t\tActive:\t\t\t%d GB'%(vmStats['Pages active']/1024/1024/1024 ))
    print('\t\tInactive:\t\t%d GB'%(vmStats['Pages inactive']/1024/1024/1024))
    print('\t\tFree:\t\t\t%d GB'%(vmStats['Pages free']/1024/1024/1024))
    print('\t\tReal Memory Total (ps):\t%.3f GB'%(rssTotal/1024/1024/1024))
    print('\t\t=====================================')
    print('---------------------------------------------------------------')


def directory_structure(study_area, data, shps, results, tests=None):
    """set up directory structure for project
    
    Parameters
    ----------
    study_area : str
        county of interest
    data : str
        upper-level data directory.
    shps : list
        subdir paths to initial.
    results : str
        upper-level results directory.
    tests : list
        names of test areas within the study area. Default is None.
    
    Returns
    -------
    ini, int, cle, cc, co, cn, cx, rpl, tbl, rfl : tuple
    """
    
    
    def _set_paths(area, data, results):
        """sets and returns strings of paths
        
        Parameters
        ----------
        area : str
            county of interest or test study area within the county.
        data : str
            upper-level data directory.
        results : str
            upper-level results directory.
        
        Returns
        -------
        ini, int, cle, cc, co, cn, cx, ca, \
        rpl, tbl, rfl, clean_ss, results_ss : tuple
        """
        
        # data directories
        d = data + area + '/'
        ini, int, cle = d + 'initial/', d + 'intermediary/', d + 'clean/'
        
        # clean data sub directories
        cc, co = cle + 'census_data/', cle + 'observation_data/'
        cn, cx = cle + 'network_data/', cle + 'cost_matrices/'
        ca = cle + 'allocation_data/'
        clean_ss = [cc, co, cn, cx, ca]
        
        # results directories
        results = results + area + '/'
        
        # results sub directories
        rpl = results + 'plots/'
        tbl = results + 'tables/'
        rfl = results + 'facility_location/'
        results_ss = [rpl, tbl, rfl]
        
        return ini, int, cle, cc, co, cn, cx, ca,\
               rpl, tbl, rfl, clean_ss, results_ss
    
    
    def _set_dirs(initi, shps, inter, clean, result):
        """create directories
        
        Parameters
        ----------
        initi : str
            path to initial data directory.
        shps : list
            subdir paths to initial.
        inter : str
            path to intermediary data directory.
        clean : list
            subdir paths to to clean data directory.
        result : list
            subdir paths to results directory.
        """
        
        for shp in shps:
            if not os.path.exists(initi+shp):
                os.makedirs(initi+shp)
        
        if not os.path.exists(inter):
            os.makedirs(inter)
        
        for sub_dir in clean:
            if not os.path.exists(sub_dir):
                os.makedirs(sub_dir)
        
        for sub_dir in result:
            if not os.path.exists(sub_dir):
                os.makedirs(sub_dir)
    
    
    # for each test case set up a directory
    if tests:
        for test in tests:
            
            # get path names
            ini, int, cle,\
            cc, co, cn, cx, ca,\
            rpl, tbl, rfl, clean_ss, results_ss = _set_paths(test,
                                                             data,
                                                             results)
            
            # create dirs
            _set_dirs(ini, shps, int, clean_ss, results_ss)
    
    
    # set up a directory for the main study area
    # get path names
    ini, int, cle,\
    cc, co, cn, cx, ca,\
    rpl, tbl, rfl, clean_ss, results_ss = _set_paths(study_area, data, results)
    
    # create dirs
    _set_dirs(ini, shps, int, clean_ss, results_ss)
    
    return ini, int, cle, cc, co, cn, cx, ca, rpl, tbl, rfl


def check_progress(num=None, denom=None, cTime=None):
    """Check the progress of a process with 10% completion
    intervals, and prints the current progress + cumulative
    time of current process.
    
    
    Parameters
    ----------
    num : int
        number of iterations complete
    denom : int
        total number of iterations
    cTime : float
        process start time
    """
    percentage_complete = str(int(round((float(num)/float(denom))*100.)))
    cumulative_time = round((time.time()-cTime)/60., 3)
    print('\t\t\t\t\t', percentage_complete+'% in', cumulative_time, 'min.')


def time_phase(phase=None, start=False, end=None, study_area=None):
    """
    Simple phase timer and progress checker
    
    Parameters
    ----------
    phase : str
        phase of the program. Default is None.
    start : bool
        initialization of the phase. Default is False.
    end : float
        end time of the phase. Default is None.
    study_area : str
        study area. Default is None.
    
    Returns
    -------
    time.time() : initial phase time
    """
    
    if start:
        if phase == '1':
            text = 'Fetch census geographies with CenPy: ' + study_area
        if phase == '2':
            text = 'Observations (Fire stations): ' + study_area
        if phase == '3':
            text = 'Create network: ' + study_area
        if phase == '4':
            text = 'Voronoi and PP2N creation: ' + study_area
        if phase == '5':
            text = 'Snapping observations to network: ' + study_area
        if phase == '6':
            text = 'Calculating cost matrices: ' + study_area
        if phase == '6.1':
            text = 'Segment-to-Segment subpahse: ' + study_area
        if phase == '7':
            text = 'Hard code plots: ' + study_area
        if phase == '8':
            text = 'Summary Stats: ' + study_area
        if phase == '9':
            text = 'Facility Location: ' + study_area
        print('Initializing Phase', phase + ':\t', text)
        return time.time()
    
    if end:
        end = round((time.time()-end)/60., 10)
        print('Phase', phase, 'complete')
        print('\tTotal minutes for phase', phase + ':\t', end)
        print('\tPhase', phase, 'complete:\t\t',\
              time.strftime('%Y-%m-%d %H:%M'))
        print('--------------------------------------------------------------')


def record_filter(df, column=None, sval=None, mval=None, oper=None):
    """used in phase 2 with incidents
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        dataframe of incident records
    oper : operator object *OR* str
        {(operator.eq, operator.ne), ('in', 'out')}
    sval : str, int, float, bool, etc.
        single value to filter
    mval : list
        multiple values to filter
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        dataframe of incident records
    """
    
    # use index or specific column
    if column == 'index':
        frame_col = df.index
    else:
        frame_col = df[column]
    
    # single value in column
    if not sval == None:
        return df[oper(frame_col, sval)].copy()
    
    # multiple values in column
    if not mval == None:
        if oper == 'in':
            return df[frame_col.isin(mval)].copy()
        if oper == 'out':
            return df[~frame_col.isin(mval)].copy()


def set_crs(df, proj_init=None, proj_trans=None, crs=None):
    """Set and transform the coordinate
    reference system of a geodataframe.
    
    Parameters
    ----------
    df : geopandas.GeoDataframe
        geodataframe being transformed
    proj_init : int
        intial coordinate reference system. default is None.
    proj_trans : int
        transformed coordinate reference system. default is None.
    crs : dict
        crs from another geodataframe
    
    Returns
    -------
    df : geopandas.GeoDataframe
        transformed geodataframe
    """
    
    if proj_init:
        df.crs = {'init': 'epsg:'+str(proj_init)}
    
    if proj_trans:
        df = df.to_crs(epsg=int(proj_trans))
    
    if crs:
        df = df.to_crs(crs)
    
    return df


def geom_to_float(df, xval=None, yval=None, geo_col=None, geom_type=None):
    """convert a geometric point object to single floats
    for inclusion in a dataframe column.
    
    Parameters
    ----------
    df : geopandas.GeoDataframe
        initial dataframe
    xval : str
        x coordinate column name. Default is None.
    yval : str
        y coordinate column name. Default is None.
    geo_col : str
        geomtery column name. Default is None.
    geom_type : str
        geometry type to transform into. Currently either cent
        (centroid) or repp (representative point).
    
    Returns
    -------
    df : geopandas.GeoDataframe
        updated dataframe
    """
    
    geoms = {'cent':'centroid', 'repp':'representative_point'} 
    
    try:
        # for centroids
        df[xval] = [getattr(p, geoms[geom_type]).x for p in df[geo_col]]
        df[yval] = [getattr(p, geoms[geom_type]).y for p in df[geo_col]]
    
    except AttributeError:
        try:
            # for representative points
            df[xval] = [getattr(p, geoms[geom_type])().x for p in df[geo_col]]
            df[yval] = [getattr(p, geoms[geom_type])().y for p in df[geo_col]]
        except:
            raise AttributeError(geoms[geom_type]+' attribute not present.')
    
    df.drop([xval, yval], axis=1, inplace=True)
    
    return df


def generate_xyid(df=None, geom_type='node', geo_col=None):
    """create a string xy id
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        geometry dataframe. Default is None.
    geom_type : str
        either node of segm. Default is 'node'.
    geo_col : str
        geomtery column name. Default is None.
    
    Returns
    -------
    xyid : list
        list of combined x-coord + y-coords strings
    """
    
    xyid = []
    
    for idx, geom in enumerate(df[geo_col]):
        
        if geom_type == 'segm':
            xys = ['x'+str(x)+'y'+str(y) for (x,y) in geom.coords[:]]
            xyid.append([idx, xys])
        
        # try to make the xyid from a polygon
        if geom_type == 'node':
            try:
                xy = 'x'+str(geom.centroid.x)+'y'+str(geom.centroid.y)
            
            # if the geometry is not polygon, but already point
            except AttributeError:
                try:
                    xy = 'x'+str(geom.x)+'y'+str(geom.y)
                except:
                    print('geom:', type(geom))
                    print(dir(geom))
                    raise AttributeError('geom has neither attribute:\n'\
                                         +'\t\t- `.centroid.[coord]`\n'\
                                         +'\t\t- `.[coord]`')
            
            xyid.append([idx,[xy]])
    
    return xyid


def fill_frame(frame, full=False, idx='index',
               col=None, data=None, add_factor=0):
    """fill a dataframe with a column of data
    
    Parameters
    ----------
    frame : geopandas.GeoDataFrame
        geometry dataframe
    full : bool
        create a new column (False) or a new frame (True).
        Default is False.
    idx : str
        index column name. Default is 'index'.
    col : str or list
         New column name(s). Default is None.
    data : list *OR* dict
        list of data to fill the column. Default is None. OR
        dict of data to fill the records. Default is None.
    add_factor : int
        used when dataframe index does not start at zero.
        Default is zero.
    
    Returns
    -------
    frame : geopandas.GeoDataFrame
        updated geometry dataframe
    """
    
    # create a full geopandas.GeoDataFrame
    if full:
        out_frame = gpd.GeoDataFrame.from_dict(data, orient='index')
        
        return out_frame
    
    # write a single column in a geopandas.GeoDataFrame
    else:
        frame[col] = np.nan
        for (k,v) in data:
            k += add_factor
            
            if col == 'CC':
                frame.loc[frame[idx].isin(v), col] = k
            
            elif idx == 'index':
                frame.loc[k, col] = str(v)
            
            else:
                frame.loc[(frame[idx] == k), col] = str(v)
        
        if col == 'CC':
            frame[col] = frame[col].astype('category').astype(int)
        
        return frame


def fs_desc_vars(fs_file, name='desc_var', file_type='.shp'):
    """hacky function to add in decision variables
    (already have firestations in RDC...)
    
    Parameters
    ----------
    fs_file : str
        full path and file name of fire stations
    name : str
        decision varibale column name. Default is `desc_var`. 
    file_type : str
        file extension. Default is '.shp'.
    """
    
    fs = gpd.read_file(fs_file+file_type)
    
    fs[name] = desc_var_list('service', fs.shape[0])
    
    fs.to_file(fs_file+file_type)


def desc_var_list(category, shape):
    """create decision variables for observation sets
    
    Parameters
    ----------
    category : str
        variable type. either client or service
    shape : int
        number of decision variable to create
    
    Returns
    -------
    dv_list : list
        all associated decision variables
    """
    
    if category == 'client':
        var = 'x'
    
    elif category == 'service':
        var = 'y'
    
    else:
        raise ValueError(category + ' not valid value for variable: category')
    
    dv_list = [var+'%i' % idx for idx in range(shape)]
    
    return dv_list


def in_census_geog(gdf, in_geog, subset_cols, proj,
                   buffer, file_type, xyid=None, desc_var=None):
    """perform overlay for extracting census geography location
    
    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        observation data.
    in_geog : list
        spatial join with census geography.
        Format is [geography column name,  path to data].
    subset_cols : list
        columns to subset following overlay.
    proj : int
        initial projection.
    buffer : {int, float}
        buffer in order to perform overlay.
    file_type : str
        file extension.
    xyid : str
        combined x-coord + y-coords string ID. Default is None.
    desc_var : str
        decision variable column name.  Default is None.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        observation data with column for `in_geog` id.
    """
    
    # read in census geographies for overlay
    geogs = gpd.read_file('%s%s' % (in_geog[1], file_type))
    geogs = set_crs(geogs, proj_init=proj)
    geogs = geogs[[geogs.geometry.name, 'GEOID', 'TRACT']]
    
    # convert geometries to micro polygons for overlay
    gdf.geometry = gdf.geometry.buffer(buffer)
    
    # intersection overlay with blocks in tract
    gdf = gpd.overlay(gdf, geogs, how='intersection')
    gdf = set_crs(gdf, proj_init=proj)
    
    rename_cols = {}
    if 'GEOID' in gdf.columns:
        rename_cols.update({'GEOID': in_geog[0]})
    if 'GEOID_1' in gdf.columns:
        rename_cols.update({'GEOID_1': 'GEOID'})
    if 'GEOID_2' in gdf.columns:
        rename_cols.update({'GEOID_2': in_geog[0]})
    
    check_col_prescence = dict(zip(['%s_1' % xyid, '%s_1' % desc_var],
                                   [xyid, desc_var]))
    for k, v in check_col_prescence.items():
        if k in gdf.columns:
            rename_cols.update({k: v})
    
    # rename GEOID column
    gdf = gdf.rename(index=str, columns=rename_cols)
    
    # convert geometries back to points
    gdf.geometry = gdf.centroid
    
    # keep only these columns
    if subset_cols:
        gdf = gdf[subset_cols]
    
    return gdf


def tract_subset(tract, full_cen_dir=None, subset_cen_dir=None,
                 full_obs_dir=None, subset_obs_dir=None, geos_2_shp=None,
                 geos_clip=None, fs_file_in=None, fs_file_out=None,
                 proj_init=None, area=None, area_file=None,
                 clip_by='Census Tracts', geom_unit='TRACT',
                 file_type='.shp'):
    """Create and write out the 'one tract subset'
    within Leon County, FL.
    
    Parameters
    ----------
    tract : str
        tract number
    full_cen_dir : str
        path to full clean census data
    subset_cen_dir : str
        path to subset clean census data
    full_obs_dir : 
        path to full clean observation data
    subset_obs_dir : 
        path to subset clean observation data
    geos_2_shp : dict
        see `runner_dissertation.py`
    geos_clip : list
        geographies to extract from with tract
    fs_file_in : str
        path to read in full clean firestation data
    fs_file_out : str
        path to save subset clean firestation data
    proj_init : int
        intial projection
    area : str
        `STUDY_AREA` from `runner_dissertation.py`
    area_file : str
        `fs_phase_file_name` from `runner_dissertation.py`
    clip_by : str
        Default is 'Census Tracts'.
    geom_unit : str
        Default is 'Tract'.
    file_type : str
        file extension. Default is '.shp'.
    """
    
    # ------------- Census geographies
    # --- polygons
    # write to here
    file_name = geos_2_shp[clip_by] + file_type
    one_tract_path = subset_cen_dir + file_name
    if not os.path.exists(one_tract_path):
        full_path = full_cen_dir + file_name
        tracts_gdf = gpd.read_file(full_path)
        one_tract_df = tracts_gdf[tracts_gdf[geom_unit] == tract]
        # set crs
        one_tract_df = set_crs(one_tract_df, proj_init=proj_init)
        # write
        one_tract_df.to_file(one_tract_path)
    else:
        one_tract_df = gpd.read_file(one_tract_path)
    
    # --- centroids
    # write to here
    centroid_path = subset_cen_dir + 'Centroid_' + file_name
    if not os.path.exists(centroid_path):
        _centroid_df = copy.deepcopy(one_tract_df)
        _centroid_df.geometry = _centroid_df.centroid
        # set crs
        _centroid_df = set_crs(_centroid_df, proj_init=proj_init)
        # write
        _centroid_df.to_file(centroid_path)
    
    # for each geometry to clip
    for cen_geog in geos_clip:
        
        if cen_geog.split('_')[0] == 'Centroid'\
        and cen_geog.split('_')[1] == 'Populated':
            
            cen_geog = '_'.join([cen_geog.split('_')[1]]\
                                 + [cen_geog.split('_')[0]]\
                                 +  cen_geog.split('_')[2:])
        
        # set geography file path
        geom_path = cen_geog + file_type
        full_path_temp = full_cen_dir + geom_path
        
        # write out file
        one_geom_path = subset_cen_dir + geom_path
        
        if not os.path.exists(one_geom_path):
            # read in .shp and subset
            geoms = gpd.read_file(full_path_temp)
            one_geom_df = geoms[geoms[geom_unit] == tract]
            # set crs
            one_geom_df = set_crs(one_geom_df, proj_init=proj_init)
            
            # write
            one_geom_df.to_file(one_geom_path)
    
    # ------------- Synthetic Fire Stations
    fs_full = gpd.read_file(fs_file_in + file_type)
    fs_place_time_list = area_file.split('_')
    fs = fs_place_time_list[0]
    place_time = '_' + '_'.join(fs_place_time_list[1:])
    gdf = create_synthetic_locs(fs_full, area, data_type='fire_stations')
    gdf = set_crs(gdf, proj_init=proj_init)
    path_to_save = '%s%s_Synthetic%s%s' % (subset_obs_dir, fs,\
                                           place_time, file_type)
    gdf.to_file(path_to_save)


def create_synthetic_locs(empir, area=None, data_type=None):
    """top-level function for creating sythetic locations
    
    Parameters
    ----------
    empir : geopandas.GeoDataFrame
        empirical location data
    area : str
        study area
    data_type : str
        either `firestations` or `households`
    
    Returns
    -------
    synth_df : geopandas.GeoDataFrame
        synthetic location data
    """
    
    Point = shapely.geometry.Point
    
    # fetch pre-processed, synthetic coordinates
    synth_locs = fetch_synth_locs(data_type=data_type, subset=area)
    
    # create shapely.Points from synthetic coordinates
    synth_points = [Point(p[0], p[1]) for p in synth_locs]
    
    # instantiate geodataframe
    synth_gdf = gpd.GeoDataFrame(geometry=synth_points)
    
    # purely synthetic data
    if isinstance(empir, type(None)):
        return synth_gdf
    
    # set coordinate reference system
    synth_gdf.crs = empir.crs
    
    # mimic empirical attribute structure in synthetic data
    synth_gdf = mimic_empirical_struct(synth_gdf, empir,
                                       data_type=data_type)
    
    return synth_gdf


def mimic_empirical_struct(synth, empir, data_type=None):
    """mid-level function for creating sythetic locations
    
    Parameters
    ----------
    synth : geopandas.GeoDataFrame
        synthetic location data
    empir : geopandas.GeoDataFrame
        empirical location data
    data_type : str
        either `firestations` or `households`
    
    Returns
    -------
    synth : geopandas.GeoDataFrame
        synthetic location data mimicing empirical
        attribute data
    """
    
    # fetch list of column names in empircal dataset
    cols = list(empir.columns)
    
    # remove original geometry column
    cols.remove(synth.geometry.name)
    
    # seperate general ids from dataset specific columns
    cols_general = ['xyid', 'desc_var']
    cols = [col for col in cols if col not in cols_general]
    
    if data_type == 'fire_stations':
        for col in cols:
            # recreate dataset specific IDs
            if col == 'STATION_NU':
                 synth[col] = ['synth_%s' % str(station+1)
                               for station in range(synth.shape[0])]
            # set as nan for all others
            else:
                synth[col] = np.nan
    
    # recreate important IDs (either client or service)
    for col in cols_general:
        if col == 'xyid':
            xyids = generate_xyid(df=synth, geom_type='node',
                                  geo_col=synth.geometry.name)
            synth[col] = [str(xy[1]) for xy in xyids]
        if col == 'desc_var':
            if data_type == 'households':
                var_type = 'client'
            else:
                var_type = 'service'
            synth[col] = desc_var_list(var_type, synth.shape[0])
    
    return synth


def get_fips(st, ct):
    """return cenus FIPS codes for states and counties
    
    Parameters
    ----------
    st : str
        state name
    ct : str
        county name
    
    Returns
    -------
    sf, cf  : tuple
        state and county fips
    """
    
    if st.lower() == 'fl' and ct.lower() == 'leon':
        sf, cf = '12', '073'
    
    else:
        fips_csv_path = '../us_county_fips_2010.csv'
        us_df = pd.read_csv(fips_csv_path, dtype={'ST_FIPS': str,
                                                  'CT_FIPS': str})
        us_df['CT'] = us_df['CT_Name'].apply(lambda x:\
                                            ''.join(x.split(' ')[:-1]).lower())
        st_df = us_df[us_df.ST_Post == st.upper()]
        record = st_df[st_df.CT == ct.replace(' ', '').lower()]
        sf, cf = record.ST_FIPS.values[0], record.CT_FIPS.values[0]
    
    return sf, cf


def discard_troublemakers(study_area, tiger_edges, tiger_roads):
    """Discard specific troublemaker line segments by at the county
    level. This is based on (1) in situ observation (ground
    truthed/local knowledge); (2) remotely-sensed observation (satelite
    imagery); (3) obvious errors in digitization. Troublemakers may also
    be [part of] road segments that exist is reality, but lead to errors
    in network instantiation.
    
    Parameters
    ----------
    study_area : str
        county and state name. Default is None.
    tiger_edges : bool
        using tiger edges file. Default is False.
    tiger_roads : bool
        using tiger roads file. Default is False.
    
    Returns
    -------
    discard_segs : list
        specific segment ids to drop.
    """
    
    if study_area == 'Leon_FL':
        
        # if using tiger roads
        if tiger_roads:
            OLD_BLAIRSTONE = ['1102179328854']
            SILVER_LAKE = ['110159184599', '110159062008', '110159184931']
            discard_segs = OLD_BLAIRSTONE + SILVER_LAKE
        
        if tiger_edges:
            OLD_BLAIRSTONE = [618799725, 618799786, 618799785, 634069404,
                              618799771, 618799763, 610197711]
            SILVER_LAKE = [82844664, 82844666, 82844213, 82844669,
                           82844657, 82844670, 82844652, 82844673]
            #TALLY_SQUARE = [634069168]
            discard_segs = OLD_BLAIRSTONE + SILVER_LAKE# + TALLY_SQUARE
            
    return discard_segs


def get_discard_mtfcc_by_desc():
    """discard these road types from the mtfcc categories
    """
    
    return ['Bike Path or Trail', 'Parking Lot Road',  'Alley',\
            'Vehicular Trail (4WD)', 'Walkway/Pedestrian Trail',\
            'Private Road for service vehicles (logging, oil fields, '\
            +'ranches, etc.)']


def fetch_synth_locs(data_type=None, subset=None):
    """mid-level function for creating sythetic locations
    
    Parameters
    ----------
    data_type : str
        either `firestations` or `households`
    subset : str
        location subset name
    
    Returns
    -------
    locs : list
        synthetic coordinates specified by location type
    """
    
    leon_fl = '_Leon_FL'
    SUBSETS = ['Test_Tract', 'Test_Grid']
    SUBSETS_MADE = [ss+'%s' % leon_fl for ss in SUBSETS]
    if not subset in SUBSETS_MADE:
        raise ValueError('Subset not created for: %s' % subset)
    
    # synthetic firestations --> these are not actual fire stations
    if data_type == 'fire_stations':
        if subset == 'Test_Tract_Leon_FL':
            locs = [[623658.0362039418, 165522.0073909023],
                     [622846.360173869, 164116.5534190332],
                     [624985.5508821681, 164956.56671589593],
                     [621425.6664735656, 165131.65812323528],
                     [621915.8577637655, 163299.8906703833],
                     [623051.0375937019, 162267.9090068047],
                     [624095.9190280752, 164925.26179051955],
                     [624637.709401454, 163506.28700309902],
                     [624586.110318275, 164254.47370919347],
                     [623205.8348432387, 164318.97256316713],
                     [622315.7506584021, 164177.07508442507],
                     [623012.3382813177, 165866.945058535],
                     [623812.1240705911, 163196.69250402544],
                     [622999.438510523, 163196.69250402544],
                     [622431.8485955547, 165183.25720641422]]
    
    # synthetic households --> these are not actual household locations
    if data_type == 'households':
        if subset == 'Test_Grid_Leon_FL':
            locs = [[6.173800738007382, 8.832702952029521],
                     [7.533763837638379, 4.561023985239852],
                     [8.789114391143915, 1.928274907749076],
                     [5.197416974169744, 7.64709409594096],
                     [3.820018450184503, 5.24100553505535],
                     [4.482564575645758, 2.695433579335792],
                     [0.22832103321033176, 8.745525830258304]]
    
    return locs


def get_mtfcc_types():
    """read in dictionary of MTFCC road type descriptions
    https://www.census.gov/geo/reference/mtfcc.html
    
    ******* Ranks are subjective *******
    
    """
   
    mtfcc_types = {'S1100':{'FClass':'Primary Road',
                             'Desc': 'Primary roads are generally divided, '\
                                     +'limited-access highways within the '\
                                     +'interstate highway system or under '\
                                     +'state management, and are '\
                                     +'distinguished by the presence of '\
                                     +'interchanges. These highways are '\
                                     +'accessible by ramps and may include '\
                                     +'some toll highways.',
                             'MTFCCRank':1,
                             'Buffer':0.},
                   'S1200':{'FClass':'Secondary Road',
                             'Desc': 'Secondary roads are main arteries, '\
                                     +'usually in the U.S. Highway, State '\
                                     +'Highway or County Highway system. '\
                                     +'These roads have one or more lanes of '\
                                     +'traffic in each direction, may or may '\
                                     +'not be divided, and usually have '\
                                     +'at-grade intersections with many '\
                                     +'other roads and driveways. They often '\
                                     +'have both a local name and a route '\
                                     +'number.',
                             'MTFCCRank':2,
                             'Buffer':100.},
                   'S1400':{'FClass':'Local Neighborhood Road, '\
                                      +'Rural Road, City Street',
                             'Desc': 'Generally a paved non-arterial street, '\
                                     +'road, or byway that usually has a '\
                                     +'single lane of traffic in each '\
                                     +'direction. Roads in this feature '\
                                     +'class may be privately or publicly '\
                                     +'maintained. Scenic park roads would '\
                                     +'be included in this feature class, '\
                                     +'as would (depending on the region of '\
                                     +'the country) some unpaved roads.',
                             'MTFCCRank':3,
                             'Buffer':50.},
                   'S1500':{'FClass':'Vehicular Trail (4WD)',
                             'Desc': 'An unpaved dirt trail where a '\
                                     +'four-wheel drive vehicle is required. '\
                                     +'These vehicular trails are found '\
                                     +'almost exclusively in very rural '\
                                     +'areas. Minor, unpaved roads usable by '\
                                     +'ordinary cars and trucks belong in '\
                                     +'the S1400 category.',
                             'MTFCCRank':10,
                             'Buffer':10.},
                   'S1630':{'FClass':'Ramp',
                             'Desc': 'A road that allows controlled access '\
                                     +'from adjacent roads onto a limited '\
                                     +'access highway, often in the form of '\
                                     +'a cloverleaf interchange. These roads '\
                                     +'are unaddressable.',
                             'MTFCCRank':4,
                             'Buffer':0.},
                   'S1640':{'FClass':'Service Drive usually along a limited '\
                                      +'access highway',
                             'Desc': 'A road, usually paralleling a limited '\
                                     +'access highway, that provides access '\
                                     +'to structures along the highway. '\
                                     +'These roads can be named and may '\
                                     +'intersect with other roads.',
                             'MTFCCRank':5,
                             'Buffer':0.},
                   'S1710':{'FClass':'Walkway/Pedestrian Trail',
                             'Desc': 'A path that is used for walking, being '\
                                     +'either too narrow for or legally '\
                                     +'restricted from vehicular traffic.',
                             'MTFCCRank':12,
                             'Buffer':0.},
                   'S1720':{'FClass':'Stairway',
                             'Desc': 'A pedestrian passageway from one level '\
                                     +'to another by a series of steps.',
                             'MTFCCRank':99,
                             'Buffer':0.},
                   'S1730':{'FClass':'Alley',
                             'Desc': 'A service road that does not generally '\
                                     +'have associated addressed structures '\
                                     +'and is usually unnamed. It is located '\
                                     +'at the rear of buildings and '\
                                     +'properties and is used for '\
                                     +'deliveries.',
                             'MTFCCRank':6,
                             'Buffer':10.},
                   'S1740':{'FClass':'Private Road for service vehicles '\
                                      + '(logging, oil fields, ranches, etc.)',
                             'Desc': 'A road within private property that is '\
                                     +'privately maintained for service, '\
                                     +'extractive, or other purposes. These '\
                                     +'roads are often unnamed.',
                             'MTFCCRank':7,
                             'Buffer':0.},
                   'S1750':{'FClass':'Internal U.S. Census Bureau use',
                             'Desc': 'Internal U.S. Census Bureau use',
                             'MTFCCRank':8,
                             'Buffer':10.},
                   'S1780':{'FClass':'Parking Lot Road',
                             'Desc': 'The main travel route for vehicles '\
                                     +'through a paved parking area.',
                             'MTFCCRank':9,
                             'Buffer':0.},
                   'S1820':{'FClass':'Bike Path or Trail',
                             'Desc': 'A path that is used for manual or '\
                                     +'small, motorized bicycles, being '\
                                     +'either too narrow for or legally '\
                                     +'restricted from vehicular traffic.',
                             'MTFCCRank':11,
                             'Buffer':0.},
                   'S1830':{'FClass':'Bridle Path',
                             'Desc': 'A path that is used for horses, being '\
                                     +'either too narrow for or legally '\
                                     +'restricted from vehicular traffic.',
                             'MTFCCRank':13,
                             'Buffer':0.},
                   'S2000':{'FClass':'Road Median',
                             'Desc': 'The unpaved area or barrier between '\
                                     +'the carriageways of a divided road.',
                             'MTFCCRank':14,
                             'Buffer':0.}
                   }
    
    return mtfcc_types
