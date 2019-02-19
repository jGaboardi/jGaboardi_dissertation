
"""

waverly blocks



"""
import numpy as np
import geopandas as gpd
import pandas as pd
import sys

from optimization_models import optimization_models as optm
print((dir(optm)))


models = ["lscp", "mclp", "pmp", "pcp"]



CO = sys.argv[1].capitalize()                   # Set county
test_cases = ["Waverly", "Indianhead", "Small", "Interstate", "Campus"]
if CO in test_cases:
    CO = "Test_"+CO
YR = "2010"                                     # Set year
GEOMETRY = "geometry"                           # set geometry column name
LENGTH = "length"                               # set length column name
NTW_SEGM_ID_NAME = "SegID"                      # ID name for network segments
NTW_NODE_ID_NAME = "NodeID"                     # ID name for network nodes
XY_ID = "xyID"                                  # xy string ID name
XVAL, YVAL = "CentX", "CentY"                   # individual x and y columns
SNAP_METHOD = ["segments"]                      # methods for snapping obsv
                                                # to network segments
SNAPPING_RESTRICTION = [1,4]                    # interstates and ramps
# data directories
data = "../data/"+CO+"/"
initial, inter, clean = data+"initial/", data+"intermediary/", data+"clean/"
# clean data sub directories
c_cen = clean+"census_data/"
c_obs = clean+"obs_data/"
c_cmx = clean+"cost_matrices/"




#Snapped_to_segments_CensusBlocksPopulatedCentroids.shp
client_file = c_cen+"Snapped_to_segments_CensusBlocksPopulatedCentroids.shp"
facility_file = c_cen+"Snapped_to_segments_CensusBlocksPopulatedCentroids.shp"


cli_df = gpd.read_file(client_file)
################
cli_df.GEOID = [str(geoid) for geoid in list(cli_df.GEOID)]
cli_df_index = list(cli_df.GEOID)
###################



fac_df = gpd.read_file(facility_file)


cli_df["cli_dv"] = ["x" + str(record+1) for record in range(cli_df.shape[0])]
fac_df["fac_dv"] = ["y" + str(record+1) for record in range(fac_df.shape[0])]

# this needs updating........................................................ string name.....
# change GEOID to str... not unicode... change early in the prep...
# may need other changes after doing this...
cli_index = [cli_df]




cli_str = "PopulatedBlockCentroids"
fac_str = "PopulatedBlockCentroids"
snap_type = "snapped_to_segments"
#PopulatedBlockCentroids_x_PopulatedBlockCentroids_snapped_to_segments_DataFrame.csv
cost_matrix_file = c_cmx+cli_str+"_x_"+fac_str+"_"+snap_type+"_DataFrame.csv"
cost_matrix_df = pd.read_csv(cost_matrix_file, header=0, index_col=0,
                             dtype={"index":str})

cost_matrix_df_cli_index = [str(idx) for idx in list(cost_matrix_df.index)]




coverage = 3.
res = optm.lscp(cli_vars=list(cli_df.cli_dv), fac_vars=list(fac_df.fac_dv),
                cij=cost_matrix, s=coverage)



'''
np.random.seed(1991)
model = "mclp"
client_count, facility_count = 7, 5
weights = np.random.randint(1, 10, (client_count,1))
# weights = np.array()

cost_matrix = np.random.uniform(0.1, 4., (client_count,facility_count))
cost_matrix = np.round(cost_matrix,  decimals=2)
# cost_matrix = np.array()

facilities = 2

# Location Set Cover Problem
if model == "lscp":
    coverage = 3.
    lp = model+"_s"+str(coverage)
    res = optm.lscp(cij=cost_matrix, s=coverage, write_lp=lp)
# Maximal Covering Location Problem
if model == "mclp":
    coverage = 1.5
    lp = model+"_p"+str(facilities)+"s"+str(coverage)
    res = optm.mclp(ai=weights, cij=cost_matrix, s=coverage,
                    p=facilities, write_lp=lp)
# p-median Problem
if model == "pmp":
    # set objective function
    lp = model+"_p"+str(facilities)
    res = optm.pmp(ai=weights, cij=cost_matrix, p=facilities, write_lp=lp)
# p-center Problem
if model == "pcp":
    # set objective function
    lp = model+"_p"+str(facilities)
    res = optm.pcp(cij=cost_matrix, p=facilities, write_lp=lp)
if model not in ["lscp", "mclp", "pcp", "pmp"]:
    raise Exception("Model name found. Check model acronym.")



print res
'''



































"""
Complete runner for the PP2N workflow
#
import os, sys, time, subprocess
from general_utils import utils
from spaghetti import spaghetti as spgh
import dissertation_workflows as dwfs

clean_census_data = False           ##################################
clean_service_data = False          #################################
create_network = True               #################################
snap_observations = True            #################################
calculate_cost_matrices = True      #############################
locate_facilities = True            ######################


# import the following only if outside the RDC
if sys.argv[1].upper() == "RDC":
    IN_RDC = True
else:
    IN_RDC = False
    if clean_census_data:
        from census_funcs import census_funcs as cf
    if clean_service_data:
        from service_funcs import service_funcs as sf
    utils.print_packages(IN_RDC)                      # List major packages


# Set defined variables and command line arguments
FLN_HARN_m = 2779                               # Florida North HARN meters
FLN_HARN_f = 2883                               # Florida North HARN feet
NAD83 = 4269                                    # North America Datum 1983
WGS84 = 4326                                    # World Geodetic System 1984
WEBMERCATOR = 3857                              # Web Mercator (Aux. Sphere)
ATTR1 = "MTFCC"                                 # Tiger Roads attribute (Prim)
ATTR1RANK = "MTFCCRank"                         # Tiger Roads attr rank (Prim)
ATTR2 = "TLID"                                  # Tiger EDGES attribute (Sec)
ST = sys.argv[2].upper()                        # Set state
CO = sys.argv[3].capitalize()                   # Set county
test_cases = ["Waverly", "Indianhead", "Small", "Interstate", "Campus"]
if CO in test_cases:
    CO = "Test_"+CO
if ATTR2 == "TLID":
    FROM_TIGER_ROADS = True
    FROM_TIGER_EDGES = False
if ATTR2 == "TLID":
    FROM_TIGER_EDGES = True
    FROM_TIGER_ROADS = False
YR = "2010"                                     # Set year
GEOMETRY = "geometry"                           # set geometry column name
LENGTH = "length"                               # set length column name
NTW_SEGM_ID_NAME = "SegID"                      # ID name for network segments
NTW_NODE_ID_NAME = "NodeID"                     # ID name for network nodes
XY_ID = "xyID"                                  # xy string ID name
XVAL, YVAL = "CentX", "CentY"                   # individual x and y columns
SNAP_METHOD = ["segments"]                      # methods for snapping obsv
                                                # to network segments
SNAPPING_RESTRICTION = [1,4]                    # interstates and ramps
# data directories
data = "../data/"+CO+"/"
initial, inter, clean = data+"initial/", data+"intermediary/", data+"clean/"
# clean data sub directories
c_cen = clean+"census_data/"
c_obs = clean+"obs_data/"
c_net = clean+"network_data/"
c_cmx = clean+"cost_matrices/"
# results directories
results = "../results/"+CO+"/"
# set up directories
utils.set_up_dirs(initial=initial, intermediary=inter, clean=clean,
                  clean_subs=[c_cen, c_obs, c_net, c_cmx],
                  results=results, )

#---------------------------------------------------------------------- Phase 1
if not IN_RDC and clean_census_data:
    phase = "1" # PHASE 1: Fetch census geographies
    phase_start = utils.time_phase(phase=phase, start=True, co=CO)
    # Set phase variables
    db, mps, pkg = "DecennialSF1", "tigerWMS_Census"+YR, "geopandas"
    geos = ["Census Tracts","Census Block Groups","Census Blocks","Counties"]
    # Call cenpy
    cf.call_cenpy(state=ST, county=CO, year=YR, database=db, map_service=mps,
                  geographies=geos, pkg=pkg, proj_init=WEBMERCATOR,
                  proj_trans=FLN_HARN_m, save_clean=c_cen, xy_id=XY_ID,
                  xval=XVAL, yval=YVAL, geo_col=GEOMETRY)
    utils.time_phase(phase=phase, end=phase_start)

#---------------------------------------------------------------------- Phase 2

if not IN_RDC and clean_service_data and CO == "Leon":
    # Fire Station cleanse ``````````````````````````````````````````````````2a
    phase = "2" # PHASE 2: Clean observation data
    phase_start = utils.time_phase(phase=phase, start=True, co=CO)
    subphase = "a"
    subphase_start = utils.time_phase(phase=phase, subphase=subphase,
                                    co=CO, start=True)
    file_name = "TFD_Stations.shp"
    file_dir = "TFDStations_2010/" + file_name
    raw_file = initial+file_dir
    phase_file = c_obs+file_name
    # clean fire stations
    sf.clean_fire_stations(in_file=raw_file, pf=phase_file,
                           proj_init=WGS84, proj_trans=FLN_HARN_m,
                           xy_id=XY_ID, xval=XVAL, yval=YVAL,
                           geo_col=GEOMETRY)
    utils.time_phase(phase=phase, subphase=subphase, end=subphase_start)
    # Fire Station subset ```````````````````````````````````````````````````2b
    subphase = "b"
    subphaseStart = utils.time_phase(phase=phase, subphase=subphase,
                                     co=CO, start=True)
    utils.create_subset(area="Test_Waverly", clean_full_cen=c_cen,
                        clean_full_obs=c_obs, full_file=phase_file,
                        area_file=file_name, proj_init=FLN_HARN_m)
    utils.create_subset(area="Test_Indianhead", clean_full_cen=c_cen,
                        clean_full_obs=c_obs, full_file=phase_file,
                        area_file=file_name, proj_init=FLN_HARN_m)
    utils.create_subset(area="Test_Small", clean_full_cen=c_cen,
                        clean_full_obs=c_obs, full_file=phase_file,
                        area_file=file_name, proj_init=FLN_HARN_m)
    utils.create_subset(area="Test_Interstate", clean_full_cen=c_cen,
                        clean_full_obs=c_obs, full_file=phase_file,
                        area_file=file_name, proj_init=FLN_HARN_m)
    utils.time_phase(phase=phase, subphase=subphase, end=subphase_start)
    # Incident cleanse ``````````````````````````````````````````````````````2c
    subphase = "c"
    subphaseStart = utils.time_phase(phase=phase, subphase=subphase,
                                     co=CO, start=True)
    file_name = "TallahasseeIncidentCalls_"+YR+".shp"
    file_dir = "TallahasseeIncidentCalls_"+YR+"/" + file_name
    raw_file = initial+file_dir
    only_pop_block = c_cen+"CensusBlocksPopulated.shp"
    phase_file = "TallahasseeResidentialIncidentCalls_"+YR+".shp"
    phase_dir = c_obs+phase_file
    # clean incident calls
    sf.clean_incid_calls(in_file=raw_file, pf=phase_dir, proj_init=FLN_HARN_f,
                         proj_trans=FLN_HARN_m, strip_cols=True, add_prop=True,
                         only_geom=True, serv_area=True, info_legend=True,
                         no_training=True, no_date=True, dup_inc_ids=True,
                         only_residential=True, no_motel=True,
                         no_bad_geocode=True, only_populated=only_pop_block,
                         save_intermediaries=inter, xy_id=XY_ID,
                         geo_col=GEOMETRY, xval=XVAL, yval=YVAL)
    utils.time_phase(phase=phase, subphase=subphase, end=phase_start)
    # Incident subset ```````````````````````````````````````````````````````2b
    subphase = "d"
    subphaseStart = utils.time_phase(phase=phase, subphase=subphase,
                                     co=CO, start=True)
    utils.create_subset(area="Test_Waverly", clean_full_cen=c_cen,
                        clean_full_obs=c_obs, full_file=phase_dir,
                        area_file=phase_file, proj_init=FLN_HARN_m)
    utils.create_subset(area="Test_Indianhead", clean_full_cen=c_cen,
                        clean_full_obs=c_obs, full_file=phase_dir,
                        area_file=phase_file, proj_init=FLN_HARN_m)
    utils.create_subset(area="Test_Small", clean_full_cen=c_cen,
                        clean_full_obs=c_obs, full_file=phase_dir,
                        area_file=phase_file, proj_init=FLN_HARN_m)
    utils.create_subset(area="Test_Interstate", clean_full_cen=c_cen,
                        clean_full_obs=c_obs, full_file=phase_dir,
                        area_file=phase_file, proj_init=FLN_HARN_m)
    utils.time_phase(phase=phase, subphase=subphase, end=subphase_start)

#---------------------------------------------------------------------- Phase 3
if create_network:
    phase = "3" # PHASE 3: Scrub TIGER roads data
    phase_start = utils.time_phase(phase=phase, start=True, co=CO)
    if ATTR2 == "TLID":
        edges_subset = {"col":"ROADFLG", "val":"Y"}
    else:
        edges_subset = None
    # prep a tiger road .shp file
    full_segms, full_nodes = inter+"NetSegms", inter+"NetNodes"
    simplified_segms = c_net+"SimplifiedSegms"
    simplified_nodes = c_net+"SimplifiedNodes"
    raw_info = {"initial": initial, "study_area": CO, "state": ST, "year": YR}
    net = spgh.SpaghettiNetwork(sdata=raw_info, inter=inter,
                                record_components=True,
                                sid_name=NTW_SEGM_ID_NAME,
                                nid_name=NTW_NODE_ID_NAME, xy_id=XY_ID,
                            
                                proj_init=NAD83, proj_trans=FLN_HARN_m,
                                proj_units="meters",
                                phase=phase,
                                largest_component=True,
                                len_col=LENGTH, geo_col=GEOMETRY,
                            tiger_edges=FROM_TIGER_EDGES,
                            tiger_roads=FROM_TIGER_ROADS,
                            edges_subset=edges_subset,
                            save_full=True,
                            full_net_segms=full_segms,
                            full_net_nodes=full_nodes,
                            simplify_network=True,
                            save_simplified=True,
                            simplified_net_segms=simplified_segms,
                            simplified_net_nodes=simplified_nodes,
                            gen_matrix=True, mtx_to_csv=c_net,
                            gen_adjmtx=True,
                            algo="dijkstra",
                            calc_stats=True)
    utils.time_phase(phase=phase, end=phase_start)

#---------------------------------------------------------------------- Phase 4
if snap_observations:
    phase = "4" # PHASE 4: Snap observation data to network
    phase_start = utils.time_phase(phase=phase, start=True, co=CO)
    segm_file = None

    try:
        type(net)
    except:
        net = spgh.load_pickled_network(c_net)
    dwfs.snap_obsvs(net, segm_file=segm_file, county=CO, phase=phase, clean=clean,
                    restrict_col=ATTR1RANK, restrict=SNAPPING_RESTRICTION,
                    snap_method=SNAP_METHOD, xval=XVAL, yval=YVAL, progress=True,
                    geo_col=GEOMETRY, xy_id=XY_ID)
    utils.time_phase(phase=phase, end=phase_start)

#---------------------------------------------------------------------- Phase 5
if calculate_cost_matrices:
    phase = "5" # PHASE 5: Calculate all cost matrices
    phase_start = utils.time_phase(phase=phase, start=True, co=CO)
    segm_file = None
    dwfs.matrix_calc(net, segm_file=segm_file, county=CO, phase=phase, clean=clean,
                     snap_method=SNAP_METHOD, progress=True,
                     df_to_csv=True)
    utils.time_phase(phase=phase, end=phase_start)
    
#---------------------------------------------------------------------- Phase 6
# Facility Location
if locate_facilities:

    dwfs.optim_loc(
"""
    
    
    
    
    
    
