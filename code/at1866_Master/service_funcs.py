"""Functions for filtering and manipulating firestation and incident
data in the PP2N worflow.
"""

# non-standard library imports
import geopandas as gpd

# project library imports
from . import utils


def clean_fire_stations(in_file=None, pf=None, proj_init=None,
                        proj_trans=None, xyid=None, xval=None, yval=None,
                        geo_col=None, print_diags=False, file_type='.shp'):
    """subset for the active, full fire station in Leon County, 2010 and
    write out the dataframe as a .shp file
    
    Parameters
    ----------
    in_file : str
        file path to raw data. Default is None.
    pf : str
        file path to cleaned data. Default is None.
    file_type : str
        file extension. Default is '.shp'.
    proj_init : int
        initial projection. Default is None.
    proj_trans  : int
        transformed projection. Default is None.
    xyid : str
        combined x-coord + y-coords string ID. Default is None.
    xval : str
        x coordinate column name. Default is None.
    yval : str
        y coordinate column name. Default is None.
    geo_col : str
        geometry column name. Default is None.
    print_diags : bool
        print observation cleaning diagnostics. Default is False.
    """
    
    gdf = gpd.read_file(in_file+file_type)
    active = ['FS%s' % idx for idx in range(1,16)]
    temp_sts = gdf[~gdf['STATION'].isin(active)].shape[0]
    actv_sts = gdf[gdf['STATION'].isin(active)].shape[0]
    
    # Keep only active stations
    if print_diags:
        print('\t**** Dropping fire stations that are not permanent ***')
        print('\t\tVolunteer and temporary fire stations:\t%s' % temp_sts)
        print('\t\tPermanent fire stations:\t\t%s' % actv_sts)
        print('\t-----------------------------------------------------')
    
    gdf = utils.record_filter(gdf, column='STATION', mval=active, oper='in')
    
    # Transform, add XY float columns, and write out
    gdf = utils.set_crs(gdf, proj_init=proj_init, proj_trans=proj_trans)
    gdf = utils.geom_to_float(gdf, geom_type='cent',
                              xval=xval, yval=yval, geo_col=geo_col)
    
    if xyid:
        gdf.reset_index(inplace=True, drop=True)
        xyid_list = utils.generate_xyid(df=gdf, geo_col=geo_col)
        gdf = utils.fill_frame(gdf, col=xyid, data=xyid_list)
    
    # drop  nonessential columns
    drop_cols = ['NAME', 'AINAME', 'TAXID', 'STATION']
    gdf.drop(drop_cols, axis=1, inplace=True)
    gdf.to_file('%s%s' % (pf, file_type))
