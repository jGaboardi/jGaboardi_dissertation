"""Functions for fetching and manipulating census
data in the PP2N worflow.
"""

# non-standard library imports
import geopandas as gpd

# project library imports
from . import utils


def fetch_census_data(query=None, year=None, database=None, pkg=None,
                      map_service=None, geographies=None, drop_columns=None,
                      save_initial=None, proj_init=None, file_type='.shp'):
    """
    (1) Set up the census api call and 
    (2) Call census api using cenpy for geographic data
    (3) write datafame out as a .shp file of
        the specified census geography
    
    Parameters
    ----------
    query : str
        SQL syntax query to query the census api for more info:
        https://github.com/ljwolf/cenpy
    year : str
        year in question. Default is None.
    database : str
        census database to be called. Default is None.
    pkg : str
        geodata package handler. currently supports geopandas.
        Default is None.
    map_service : str
        census map service to be called. Default is None.
    geographies : dict
        census geographies requested as keys with census geography
        shp names as values. Default is None.
    proj_init : int
        intial coordinate reference system. default is None.
    drop_columns : list
        trouble maker columns to drop. Default is None.
    save_initial : str
        full path name to initial file. Default is None.
    file_type  : str
        file extension. Default is '.shp'.
    """
    
    # only when outside the RDC
    import cenpy as cen 
    
    # -- set database
    # set specific database
    api_databases = {year:database+year}
    # set census year
    db = api_databases[year]
    # api connection instance
    api_conn = cen.base.Connection(db)
    # set the map service
    api_conn.set_mapservice(map_service)
    # set map service layers
    layers = list(api_conn.mapservice.layers.items())
    
    # key code for ESRI layer with layer name
    geokey_2_layername = {k:str(v._name) for k,v in layers
                          if list(geographies.keys()).__contains__(v._name)}
    
    # key code for ESRI layer with shapefile name to save out
    geokey_2_shpname = {k:geographies[v] for k,v\
                                         in list(geokey_2_layername.items())}
    
    for geo_key, shp_name in list(geokey_2_shpname.items()):
        
        # send query and retrieve dataframe
        df = api_conn.mapservice.query(layer=geo_key, where=query, pkg=pkg)
        
        if drop_columns:
            if set(drop_columns).intersection(set(df.columns)):
                df.drop(drop_columns, axis=1, inplace=True)
        
        df = utils.set_crs(df, proj_init=proj_init)
        
        if save_initial:
            file_out = '%s%s/%s' % (save_initial, shp_name, shp_name)
            write_file(df, df_geo=file_out, file_type=file_type)
        
        if not save_initial:
            
            return df


def prepare_census_data(geo=None, year=None, read_initial=None, proj_init=None,
                        proj_trans=None, xval=None, yval=None, geo_col=None,
                        xyid=None, save_clean=None, transform=None,
                        transform_shp=None, col_as_type=None, subsets=None,
                        desc_var=None, file_type='.shp'):
    """project .shp, add/clean columns, and create
    centroid files, write.
    
    Parameters
    ----------
    geo : str
        census geography to be prepped. Default is None.
    year : str
        year in question. Default is None.
    read_initial  : str
        full path name to initial file. Default is None.
    file_type : str
        file extension. Default is '.shp'.
    proj_init : int
        intial coordinate reference system. default is None.
    proj_trans : str
        transformed coordinate reference system. Default is None.
    save_clean : str
        full path name to clean file. Default is None.
    xyid : str
        combined x-coord + y-coords string ID. Default is None.
    desc_var : str
        decision variable column name. Default is None.
    xval : str
        x coordinate column name. Default is None.
    yval : str
        y coordinate column name. Default is None.
    geo_col : str
        geometry column name. Default is None.
    transform : str
        geometric transformation. [Option] 'centroid'. Default is None.
    transform_shp : str
        geometric transformation .shp name.  Default is None.
    col_as_type : dict
        {column: type} 
    subsets : list
        list of [name, {'relate':operator,
                        'col':pandas.DataFrame column name
                        'val': rhs of boolean indexer}]
    """
    
    df = gpd.read_file('%s%s/%s%s' % (read_initial, geo, geo, file_type))
    
    # change dataframe columns to a specified type
    if col_as_type:
        for c,t in list(col_as_type.items()):
            df[c] = df[c].astype(t)
    df = utils.set_crs(df, proj_init=proj_init, proj_trans=proj_trans)
    
    # Add centroid and rep point XY float columns to the dataframe
    df = utils.geom_to_float(df, geom_type='cent',
                             xval=xval, yval=yval, geo_col=geo_col)
    
    # add desc var
    df[desc_var] = utils.desc_var_list('client', df.shape[0])
    
    # Write polygon files
    write_file(df.copy(), df_geo='%s%s' % (save_clean, geo),
               xyid=xyid, geo_col=geo_col, file_type=file_type)
    
    if transform:
        # Write centroid files
        cent_df = df.copy()
        cent_df.geometry = cent_df.centroid
        write_file(cent_df, df_geo='%s%s' % (save_clean, transform_shp),
                   xyid=xyid, geo_col=geo_col, file_type=file_type)
    
    # Data subsets
    if subsets:
        for ss in subsets:
            df = df[ss[1]['relate'](df[ss[1]['col']], ss[1]['val'])].copy()
        
        # name
        geo = '%s_%s' %  (ss[0], geo)
        
        # rewrite desc var for subset
        df['desc_var'] = utils.desc_var_list('client', df.shape[0])
        
        # Write polygon files
        write_file(df.copy(), df_geo='%s%s' % (save_clean, geo),
                   file_type=file_type, xyid=xyid, geo_col=geo_col)
        
        if transform:
            # Write centroid files
            cent_df = df.copy()
            cent_df.geometry = cent_df.centroid
            write_file(cent_df, xyid=xyid, geo_col=geo_col,
                       df_geo='%s%s_%s' % (save_clean,ss[0], transform_shp),
                       file_type=file_type)


def write_file(df, df_geo=None, xyid=None, geo_col=None, file_type=None):
    """write out geopandas.GeoDataFrame as a shapefile
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        geometry dataframe
    df_geo : str
        file path. Default is None.
    xyid : str
        combined x-coord + y-coords string ID. Default is None.
    geo_col : str
        geometry column name. Default is None.
    file_type : str
        file extension. Default is None.
    """
    
    if xyid:
        df.reset_index(inplace=True, drop=True)
        xyid_list = utils.generate_xyid(df=df, geo_col=geo_col)
        df = utils.fill_frame(df, col=xyid, data=xyid_list)
    
    df.to_file('%s%s' % (df_geo.replace(' ', ''), file_type))


def clean_parcels(in_file=None, pf=None, proj_init=None, proj_trans=None,
                  geo_col=None, xyid=None, desc_var=None, print_diags=True,
                  in_geogs=None, min_thresh=0.0, buffer=0.001,
                  popcol=None, subset_cols=None, file_type='.shp'):
    """subset populated block / parcels overlay
    
    Parameters
    ----------
    in_file : str
        file path to raw data. Default is None.
    pf : str
        file path to cleaned data. Default is None.
    proj_init : int
        initial projection. Default is None.
    proj_trans  : int
        transformed projection. Default is None.
    geo_col : str
        geometry column name. Default is None.
    xyid : str
        combined x-coord + y-coords string ID. Default is None.
    in_geogs : list
        spatial join with census geography. Default is None.
        Format is [geography column name,  path to data].
    desc_var : str
        decision variable column name. Default is None.
    print_diags : bool
        print observation cleaning diagnostics. Default is False.
    min_thresh : {int, float}
        minimum threshold for being considered populated.
        Default is 0.0.
    buffer : {int, float}
        buffer in order to perform overlay. Default is 0.001.
    pop_col : str
        summed population estimate column. Default is None.
    subset_cols : list
        columns to subset following overlay. Default is None.
    file_type : str
        file extension. Default is '.shp'.
    """
    
    # read in all parcels
    df = gpd.read_file('%s%s' % (in_file, file_type))
    
    df = df[subset_cols + [df.geometry.name]]
    
    # transform crs
    df = utils.set_crs(df, proj_init=proj_init, proj_trans=proj_trans)
    
    # subset to only populated parcels
    pop_df = df[df[popcol] > min_thresh]
    
    # Keep only populated parcels
    if print_diags:
        print('\t**** Dropping unpopulated parcels ***')
        print('\t\tTotal parcels:\t\t\t%s' % df.shape[0])
        print('\t\tPopulated parcels:\t\t%s' % pop_df.shape[0])
        print('\t-----------------------------------------------------')
    
    # derive centroids
    gdf = pop_df.copy()
    gdf[geo_col] = gdf.centroid
    
    # add desc var
    gdf[desc_var] = utils.desc_var_list('client', gdf.shape[0])
    
    if xyid:
        gdf.reset_index(inplace=True, drop=True)
        xyid_list = utils.generate_xyid(df=gdf, geo_col=geo_col)
        gdf = utils.fill_frame(gdf, col=xyid, data=xyid_list)
    
    # perform overlay for extracting census geography location
    for in_geog in in_geogs:
        gdf = utils.in_census_geog(gdf, in_geog, None, proj_trans,
                                   buffer, file_type, xyid=xyid,
                                   desc_var=desc_var)
    
    # write
    gdf.to_file('%s%s' % (pf, file_type))


