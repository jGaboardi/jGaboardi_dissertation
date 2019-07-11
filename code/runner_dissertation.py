'''
Complete runner for dissertation workflow
'''

# standard library imports
import os, sys, time, subprocess, operator

# project library imports
from at1866_Master import utils
from at1866_Master import spaghetti as spgh
from at1866_Master import census_funcs as cf
from at1866_Master import dissertation_workflows as dwfs
from at1866_Master import figure_plotter as figplot
from at1866_Master import summary_stats as sumstat

print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
print('::::\tRun initialized at:\t'+time.strftime('%Y-%m-%d %H:%M')+'\t::::')
import cpuinfo
print(cpuinfo.get_cpu_info()['brand'])
print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')

fetch_census_data = False
clean_census_data = False
clean_service_data = True

# import the following only if outside the RDC
if sys.argv[1].upper() == 'RDC':
    IN_RDC = True
else:
    IN_RDC = False
    if clean_service_data:
        from at1866_Master import service_funcs as sf

# Set defined variables and command line arguments
FLN_HARN_m = 2779                           # Florida North HARN meters
NAD83 = 4269                                # North America Datum 1983
WGS84 = 4326                                # World Geodetic System 1984
WEBMERCATOR = 3857                          # Web Mercator (Aux. Sphere)
ATTR1 = 'MTFCC'                                 # Tiger attribute (Prim)
ATTR1RANK = 'MTFCCRank'                         # Tiger attr rank (Prim)
ATTR2 = 'TLID'                                  # Tiger EDGES attribute (Sec)
if ATTR2 == 'LINEARID':                         # Tiger ROADS attribute (Sec)
    TIGER_ROADS = True
    TIGER_EDGES = False
if ATTR2 == 'TLID':
    TIGER_EDGES = True
    TIGER_ROADS = False

ST = sys.argv[2].upper()                        # Set state
CO = sys.argv[3].capitalize()                   # Set county
YR = '2010'                                     # Set year
place_time = '_%s_%s_%s' % (CO, ST, YR)

 # parallel offset for pp2n
try:
    PP2N_OFFSET = float(sys.argv[4])
except IndexError:
    PP2N_OFFSET = 5.
LVOR_OFFSET = 1.                        # symdiff offset for va2n
LVOR_PNT_RHO = 300                      # LVD point density for va2n

# when running a subset
try:
    SS = 'Test_%s' % sys.argv[5].capitalize()    # Set subset within Leon County
    STUDY_AREA = '%s_%s_%s' % (SS, CO, ST)
    test_cases, test_subsets = None, None
# when running full county
except IndexError:
    SS = None
    STUDY_AREA = '%s_%s' % (CO, ST)
    if STUDY_AREA == 'Leon_FL':
        test_cases = ['Tract', 'Grid']
        test_subsets = ['Test_%s_%s' % (case, STUDY_AREA)\
                                    for case in test_cases]
    else:
        test_cases, test_subsets = None, None

# column names
GEOMETRY = 'geometry'                   # set geometry column name
LENGTH = 'length'                       # set length column name
NTW_SEGM_ID_NAME = 'SegID'              # ID name for network segments
NTW_NODE_ID_NAME = 'NodeID'             # ID name for network nodes
XYID = 'xyid'                           # xy string ID name
DV = 'desc_var'                         # decision variable ID
XVAL, YVAL = 'CentX', 'CentY'           # individual x and y columns
# segment welding and splitting
INTRST = 'S1100'                        # interstates mtfcc code
RAMP =  'S1630'                         # ramp mtfcc code
SERV_DR = 'S1640'                       # service drive mtfcc code
SPLIT_GRP = 'FULLNAME'                  # grouped by this variable
SPLIT_BY = [RAMP, SERV_DR]              # split interstates by ramps & service
SKIP_RESTR = True                       # no weld retry if still MLS
# snapping
SNAP_METHOD = 'segments'                # obsv snapping method to segments
REPRESENTATION = ['pc2n',               # traditional centroids to network
                  'pp2n',               # pp2n method
                  'va2n']               # LVD/Overlay Okabe allocation
SNAP_RESTRICT = [INTRST, RAMP, SERV_DR] # interstates, ramps, service roads

# .shp name to use to for each stage
geos = ['Counties', 'Census Tracts', 'Census Block Groups', 'Census Blocks']
shp_names = [g.replace(' ', '') for g in geos]
shp_names = ['%s_%s_%s_%s' % (g, CO, ST, YR) for g in shp_names]
geos_2_shp = dict(list(zip(geos, shp_names)))
poly_names = shp_names + ['Populated_'+shp_names[-1]]
cent_names = ['Centroid_%s' % geom_name for geom_name in poly_names]
va2n_names = ['va2n_%s' % geom_name for geom_name in poly_names]
pp2n_names = ['pp2n_%s' % geom_name for geom_name in poly_names]

# set up directory structure
data_dir = '../data/'
results_dir = '../results/'
initial, inter, clean,\
c_cen, c_obs, c_net, c_cmx,\
c_alc, r_plots, r_tables,\
r_facloc, = utils.directory_structure(STUDY_AREA, data_dir, shp_names,
                                      results_dir, tests=test_subsets)

#--------------------------------------------------------------- Phase 1
phase = '1' # PHASE 1: Clean census data
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)
#'''
# run this code block to fetch PUBLIC census data
if not IN_RDC and fetch_census_data and not SS:
    # FIPS dictionary
    state_fips, county_fips = utils.get_fips(ST, CO)
    # SQL-style query for cenpy
    query = 'STATE=%s and COUNTY=%s' % (state_fips, county_fips)
    # set database, mapping service, and geo-data package
    db, mps, pkg = 'DecennialSF1', 'tigerWMS_Census'+YR, 'geopandas'
    # Undefined columns to drop
    drop_columns = ['ACT', 'SUFFIX', 'TABBLKSUFX2', 'VINTAGE']
    # Call cenpy
    cf.fetch_census_data(query=query, year=YR, database=db, pkg=pkg,
                         map_service=mps, geographies=geos_2_shp,
                         proj_init=WEBMERCATOR, drop_columns=drop_columns,
                         save_initial=initial)
#'''
########################################################################
# run this code block to fetch RESTRICTED census data
#add in function call for RDC data prep here....
########################################################################

#'''
# run this code block to clean PUBLIC census data
if clean_census_data and not SS:
    # change datatype of these columns
    col_as_type = {'GEOID': str}
    # prepare_census for each geography
    for idx, shp in enumerate(shp_names):
        # name of the geometric transformation .shp file
        transform_shp = cent_names[idx]
        if shp.split('_')[0] == 'CensusBlocks':
            subsets = [['Populated', {'relate': operator.gt,
                                      'col':'POP100', 'val':0}]]
        else:
            subsets = None
        # transform and subset census data
        cf.prepare_census_data(geo=shp, read_initial=initial, year=YR,
                               proj_init=WEBMERCATOR, proj_trans=FLN_HARN_m,
                               xval=XVAL, yval=YVAL, geo_col=GEOMETRY,
                               xyid=XYID, save_clean=c_cen, subsets=subsets,
                               transform='Centroid', col_as_type=col_as_type,
                               transform_shp=transform_shp, desc_var=DV)
    
#'''
utils.time_phase(phase=phase, end=phase_start)

#--------------------------------------------------------------- Phase 2
########################################################################
# Cut 'phase 2: cleaning out once ready to migrate code base############
# Keep phase 2: subset #################################################
########################################################################
phase = '2' # PHASE 2: Clean observation data
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)
cen_dir = '/clean/census_data/'
obs_dir = '/clean/observation_data/'

# fire stations file paths
fs_file_base = 'FireStations'
fs_file_name = '%s_%s_%s' % (fs_file_base, CO, ST)
fs_file_dir = '%s/%s' % (fs_file_name, fs_file_name)
fs_raw_file = '%s%s' % (initial, fs_file_dir)
fs_phase_file_name = '%s%s' % (fs_file_base, place_time)
fs_phase_file = '%s%s' % (c_obs, fs_phase_file_name)

# weighted parcel file
parcels_name = 'WeightedParcels'
parcel_file_base = '%s%s' % (parcels_name, place_time)
parcel_file_dir = '%s/%s' % (parcel_file_base, parcel_file_base)
parcel_raw_file = '%s%s' % (initial, parcel_file_dir)
parcel_phase_file = '%s%s' % (c_cen, parcel_file_base)

#########################################################  BELOW non-RDC
#'''
if not IN_RDC and clean_service_data and STUDY_AREA == 'Leon_FL':
    # Fire Station cleanse
    sf.clean_fire_stations(in_file=fs_raw_file, pf=fs_phase_file,
                           proj_init=WGS84, proj_trans=FLN_HARN_m,
                           xyid=XYID, xval=XVAL, yval=YVAL,
                           geo_col=GEOMETRY, print_diags=True)
    
    parpopcol = 'SUM_EST_PO'
    subset_cols = ['PARCEL_ID', parpopcol]
    in_geogs = [['in_block', '%s%s' % (c_cen, shp_names[-1])]]

    # Weighted Parcel
    parcels_exist = os.path.exists('%s%s' % (parcel_phase_file, '.shp'))
    if not parcels_exist:
        cf.clean_parcels(in_file=parcel_raw_file, pf=parcel_phase_file,
                         proj_init=WGS84, proj_trans=FLN_HARN_m,
                         geo_col=GEOMETRY, xyid=XYID, desc_var=DV,
                         popcol=parpopcol, subset_cols=subset_cols,
                         in_geogs=in_geogs,min_thresh=0.0)

#'''
#########################################################  ABOVE non-RDC
#'''
# hacky step to add in decision variables
# (already have firestations in RDC...)
if clean_service_data and STUDY_AREA == 'Leon_FL':
    utils.fs_desc_vars(fs_phase_file, name=DV)
#'''

# Make Tract subset
if STUDY_AREA == 'Test_Tract_Leon_FL':
    # tract of interest
    toi = '001700'
    # reset census directory and file name
    full_area_file = '_'.join(place_time.split('_')[1:3])
    full_cen_dir = '%s%s%s' % (data_dir, full_area_file, cen_dir)
    full_obs_dir = '%s%s%s' % (data_dir, full_area_file, obs_dir)
    # reset full input observations file names
    __fs_phase_file = '%s%s' % (full_obs_dir, fs_phase_file_name)
    # extract these geographies from with 'tract'
    geos_clip = poly_names[2:] + cent_names[1:] + [parcel_file_base]
    geos_clip += ['ParcelPolygons_Leon_FL_2010']
    
    # create one tract census geographies subset
    utils.tract_subset(toi, full_cen_dir=full_cen_dir, subset_cen_dir=c_cen,
                       full_obs_dir=full_obs_dir, subset_obs_dir=c_obs,
                       geos_2_shp=geos_2_shp, geos_clip=geos_clip,
                       fs_file_in=__fs_phase_file, fs_file_out=fs_phase_file,
                       area_file=fs_phase_file_name,
                       proj_init=FLN_HARN_m, area=STUDY_AREA)

# Make Grid subset
if STUDY_AREA == 'Test_Grid_Leon_FL':
    # generate and write out synthetic grid-based data
    dwfs.generate_grid(STUDY_AREA, place_time, initial, c_cen, xyid=XYID,
                       geo_col=GEOMETRY, sid_name=NTW_SEGM_ID_NAME,
                       proj_init=WGS84)

# Make Sine subset
if STUDY_AREA == 'Test_Sine_Leon_FL':
    # generate and write out synthetic sine segment data
    dwfs.generate_sine_lines(STUDY_AREA, place_time, initial,
                             geo_col=GEOMETRY, sid_name=NTW_SEGM_ID_NAME,
                             proj_init=WGS84)

utils.time_phase(phase=phase, end=phase_start)

#--------------------------------------------------------------- Phase 3
#'''
phase = '3' # PHASE 3: Scrub TIGER roads data
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)
try:
    net = spgh.load_pickled(c_net, pickle_name='Network')
    
except FileNotFoundError:
    if ATTR2 == 'TLID':
        edge_subsets = [['edge', {'relate': operator.eq,
                                   'col':'ROADFLG', 'val':'Y'}]]
    else:
        edges_subset = None
    # prep a tiger road .shp file
    full_segms = inter + 'NetSegms' + place_time
    full_nodes = inter + 'NetNodes' + place_time
    simplified_segms = c_net+'SimplifiedSegms' + place_time
    simplified_nodes = c_net+'SimplifiedNodes' + place_time
    # MTFCC types and descriptions dictionary
    mtfcc_types = utils.get_mtfcc_types()
    # Remove these types of roads from dataset (by their description)
    mtfcc_discard = utils.get_discard_mtfcc_by_desc()
    # manual segment removal
    discard_segs = None
    if not SS:
        discard_segs = utils.discard_troublemakers(STUDY_AREA,
                                                   TIGER_EDGES, TIGER_ROADS)
        proj_init = NAD83
        proj_trans = FLN_HARN_m
        
    if STUDY_AREA == 'Test_Tract_Leon_FL':
        proj_init = FLN_HARN_m
        proj_trans = FLN_HARN_m
    
    if STUDY_AREA in ['Test_Grid_Leon_FL', 'Test_Sine_Leon_FL']:
        proj_init = WGS84
        proj_trans = WGS84
    
    if STUDY_AREA == 'Test_Sine_Leon_FL':
        largest_component = False
        simplify = False
        gen_adjmtx = False
        gen_matrix = False
        save_simplified = False
        calc_stats = False
        remove_gdfs = False
    else:
        largest_component = True
        simplify = True
        gen_adjmtx = False
        gen_matrix = True
        save_simplified = True
        calc_stats = True
        remove_gdfs = True
    
    # create network
    net = spgh.SpaghettiNetwork(segmdata=initial, sid_name=NTW_SEGM_ID_NAME,
                                nid_name=NTW_NODE_ID_NAME, proj_init=proj_init,
                                proj_trans=proj_trans, proj_units='meters',
                                inter=inter, attr1=ATTR1, attr1rank=ATTR1RANK,
                                attr2=ATTR2, study_area=STUDY_AREA, county=CO,
                                state=ST, year=YR, place_time=place_time,
                                mtfcc_types=mtfcc_types, geo_col=GEOMETRY,
                                mtfcc_discard=mtfcc_discard, xyid=XYID,
                                discard_segs=discard_segs, len_col=LENGTH,
                                tiger_edges=TIGER_EDGES, mtfcc_split=INTRST,
                                edge_subsets=edge_subsets, mtfcc_intrst=INTRST,
                                mtfcc_split_grp=SPLIT_GRP, mtfcc_ramp=RAMP,
                                mtfcc_split_by=SPLIT_BY, mtfcc_serv=SERV_DR,
                                skip_restr=SKIP_RESTR, tiger_roads=TIGER_ROADS,
                                record_components=True, record_geom=True,
                                remove_gdfs=remove_gdfs, calc_len=True, 
                                largest_component=largest_component, 
                                save_full=True,
                                full_net_segms=full_segms,
                                full_net_nodes=full_nodes, 
                                calc_stats=calc_stats,
                                simplify=simplify,
                                save_simplified=save_simplified,
                                simp_net_segms=simplified_segms,
                                simp_net_nodes=simplified_nodes,
                                gen_matrix=gen_matrix, mtx_to_csv=c_net,
                                gen_adjmtx=gen_adjmtx, algo='dijkstra')
    spgh.dump_pickled(net, c_net, pickle_name='Network')
utils.time_phase(phase=phase, end=phase_start)

#'''

#--------------------------------------------------------------- Phase 4
phase = '4' # PHASE 4: Create Voronoi and PP2N polygons
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)

# -- set geographic and non-geographic data -------------
# all geographic units to run
if STUDY_AREA == 'Test_Grid_Leon_FL':
    small_unit = 'synth_households'
elif STUDY_AREA in ['Test_Tract_Leon_FL', 'Leon_FL'] and not IN_RDC:
    small_unit = 'parcels'
elif STUDY_AREA in ['Test_Tract_Leon_FL', 'Leon_FL'] and IN_RDC:
    small_unit = 'households'
elif STUDY_AREA == 'Test_Sine_Leon_FL':
    small_unit = None
else:
    raise RunTimeError('`small_unit` not defined.')

geographic_units = ['tracts', 'block_groups',
                    'pop_blocks', small_unit]
non_geographic_units = ['FireStations']

#
proj_init = FLN_HARN_m

if STUDY_AREA == 'Test_Tract_Leon_FL':
    non_geographic_units[0] = non_geographic_units[0] + 'Synthetic'

if STUDY_AREA == 'Test_Grid_Leon_FL':
    PP2N_OFFSET = .25
    LVOR_OFFSET = .01
    LVOR_PNT_RHO = 30
    non_geographic_units = []
    geographic_units = geographic_units[-2:]
    proj_init = WGS84

if STUDY_AREA == 'Test_Sine_Leon_FL':
    LVOR_OFFSET = .01
    LVOR_PNT_RHO = 3000
    proj_init = WGS84
    
    dwfs.sine_voronoi(net, alloc_dir=c_alc,
                      vor_offset=LVOR_OFFSET, vor_rho=LVOR_PNT_RHO)
    
    figplot.ch2_lvd_pp2n(LVOR_PNT_RHO)


# -- set high-precision representation methods
representation_methods = REPRESENTATION[1:]
unit_representation = geographic_units[:-1]

# create high precision network allocation models
dwfs.allocate(net, representation_methods, segm_file=None,
              pp2n_offset=PP2N_OFFSET, vor_offset=LVOR_OFFSET,
              clean=clean, net_dir=c_net, vor_rho=LVOR_PNT_RHO,
              cen_dir=c_cen, alloc_dir=c_alc, inter=inter,
              geographic_units=unit_representation, xyid=XYID,
              restrict=SNAP_RESTRICT, geo_col=GEOMETRY,
              proj_init=proj_init, sid_name=NTW_SEGM_ID_NAME,
              restrict_col='SegMTFCC', desc_var='desc_var')

utils.time_phase(phase=phase, end=phase_start)
#'''

#--------------------------------------------------------------- Phase 5
phase = '5' # PHASE 5: Snap data to network
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)

segm_file = None
try:
    net
except NameError:
    net = spgh.load_pickled_network(c_net)

# -- set high-precision representation methods
unit_representation = geographic_units

dwfs.snap_obsvs(net, segm_file=segm_file, clean=clean,
                geographic_units=unit_representation,
                non_geographic_units=non_geographic_units,
                restrict_col=ATTR1, restrict=SNAP_RESTRICT,
                snap_to=SNAP_METHOD, representation=REPRESENTATION,
                pp2n=PP2N_OFFSET, va2n=LVOR_PNT_RHO, small_unit=small_unit)
utils.time_phase(phase=phase, end=phase_start)
#'''

#'''
#--------------------------------------------------------------- Phase 6
phase = '6' # PHASE 6: Calculate all cost matrices
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)

segm_file = None
try:
    net
except NameError:
    net = spgh.load_pickled_network(c_net)

dwfs.matrix_calc(net, segm_file=segm_file, clean=clean, df_to_csv=True,
                 mtx_to_csv=False, nearest_to_pickle=False, in_rdc=IN_RDC,
                 snap_to=SNAP_METHOD, representation=REPRESENTATION,
                 geographic_units=unit_representation, small_unit=small_unit,
                 non_geographic_units=non_geographic_units,
                 pp2n=PP2N_OFFSET, va2n=LVOR_PNT_RHO)
utils.time_phase(phase=phase, end=phase_start)
#'''

#--------------------------------------------------------------- Phase 6.1
phase = '6.1' # PHASE 6.1: segment to segment
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)

if STUDY_AREA == 'Leon_FL':
    sa_ = None
    to_fs = True
else: 
    sa_ = STUDY_AREA
    #to_fs = True############################################################################
    to_fs = False##########################################################################
    
dwfs.segment_midpoints(area_prefix=sa_, obs_dir=c_obs, cen_dir=c_cen,
                       net_dir=c_net, mtx_dir=c_cmx, alc_dir=c_alc,
                       restrict=SNAP_RESTRICT, mtfcc=ATTR1,
                       sid=NTW_SEGM_ID_NAME, xyid=XYID, geo_col=GEOMETRY,
                       to_fs=to_fs)
utils.time_phase(phase=phase, end=phase_start)

#'''
#--------------------------------------------------------------- Phase 7
phase = '7' # PHASE 7: Hard code plots
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)

if STUDY_AREA == 'Test_Grid_Leon_FL':
    figplot.ch2_grid_plots()

if STUDY_AREA == 'Test_Tract_Leon_FL':
    figplot.ch2_tract_plots(to_fs=to_fs)

if STUDY_AREA == 'Leon_FL':
    figplot.ch3_plots(to_fs=to_fs)
    figplot.ch4_plots(to_fs=to_fs)
    

utils.time_phase(phase=phase, end=phase_start)
#'''

#'''
#--------------------------------------------------------------- Phase 8
phase = '8' # PHASE 8: Summary Stats
phase_start = utils.time_phase(phase=phase, start=True, study_area=STUDY_AREA)

if STUDY_AREA == 'Test_Grid_Leon_FL':
    sumstat.ch2_table1()

if STUDY_AREA == 'Test_Tract_Leon_FL':
    sumstat.ch2_table2_table3(STUDY_AREA, chapter='2')

if STUDY_AREA == 'Leon_FL':
    sumstat.ch2_table2_table3(STUDY_AREA, chapter='3')

utils.time_phase(phase=phase, end=phase_start)


print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
print(('::::\tRun finalized at:\t'+time.strftime('%Y-%m-%d %H:%M')+'\t::::'))
print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
