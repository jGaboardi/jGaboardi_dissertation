"""pysal.spaghetti
Building network objects from TIGER/Line Edge data
"""

__author__ = 'James D. Gaboardi <jgaboardi@gmail.com>'

# standard library imports
import os, pickle, platform, re, time
from ast import literal_eval

# non-standard library imports
import geopandas as gpd
import pandas as pd
import numpy as np
from libpysal import cg
from scipy.spatial import distance_matrix

# project library imports
from . import utils
from . import sauce


class SpaghettiNetwork:
    """Instantiate a spaghetti.Spaghetti network object
    """
    
    def __init__(self, file_type='.shp', tnid='TNID', tnidf='TNIDF',
                 tnidt='TNIDT', network_instance=None, segmdata=None,
                 nodedata=None, sid_name=None, nid_name=None,
                 proj_init=None, proj_trans=None, proj_units=None, inter=None,
                 attr1=None, attr1rank=None, attr2=None, study_area=None,
                 county=None, state=None, year=None, place_time=None,
                 mtfcc_types=None, mtfcc_discard=None, discard_segs=None,
                 xyid=None, geo_col=None, len_col=None, tiger_edges=False,
                 edge_subsets=None, mtfcc_split=None, mtfcc_intrst=None,
                 mtfcc_ramp=None, mtfcc_serv=None, mtfcc_split_by=None,
                 mtfcc_split_grp=None, skip_restr=False, tiger_roads=False,
                 calc_len=False, record_components=False, record_geom=False,
                 largest_component=False, calc_stats=False, gen_matrix=False,
                 mtx_to_csv=None, gen_paths=False, paths_to_csv=None,
                 gen_adjmtx=False, adjmtx_to_csv=None, algo=None,
                 def_graph_elems=False, simplify=False, save_full=False,
                 full_net_segms=None, full_net_nodes=None,
                 save_simplified=False, simp_net_segms=None,
                 simp_net_nodes=None, remove_gdfs=False):
        
        """
        Parameters
        ----------
        file_type : str
            file extension. Default is '.shp'.
        tnid : str
            TIGER/Line node ID variable used for working with
            TIGER/Line edges. Default is 'TNID'.
        tnidf : str
            TIGER/Line 'From Node' variable used for building topology
            in TIGER/Line edges. Default is 'TNIDF'.
        tnidt : str
             TIGER/Line 'To Node' variable used for building topology in
             TIGER/Line edges. Default is 'TNIDT'.
        segmdata : str **OR** geopandas.GeoDataFrame
            path to segments data or a dataframe itself.
        nodedata : str **OR** geopandas.GeoDataFrame
            nodes data. Default is None.
        sid_name : str
            segment column name. Default is None.
        nid_name : str
            node column name. Default is None.
        proj_init : int
            initial projection. Default is None.
        proj_trans : int
            transformed projection. Default is None.
        proj_units : str
            unit of transformed projection. Default is None.
        attr1 : str
            auxillary variable being used. Default is None.
        att1rank : int
            auxillary variable rank. Default is None.
        attr2 : str
            auxillary variable being used. Either 'TLID' for tiger edges
            or 'LINEARID' for tiger roads. Default is None.
        inter : str
            file path to intermediary data. Default is None.
        study_area : str
            study area within county. Default is None.
        county : str
            county of interest. Default is None.
        state : str
            state of interest. Default is None.
        year : str
            data collection year. Default is None.
        place_time : str
            place and time descriptor. Default is None. e.g.
            '_Leon_FL_2010'
        mtfcc_types : dict
            MTFCC road type descriptions. Default is None.
            from [utils.get_mtfcc_types()]
        mtfcc_discard : list
            MTFCC types (by code) to discard. Default is None.
            from [utils.get_discard_mtfcc_by_desc()]
        discard_segs : list
            specifc segment ids to discard. Default is None.
            from [utils.discard_troublemakers()]
        xyid : str
            combined x-coord + y-coords string ID. Default is None.
        geo_col : str
            geometry column name. Default is None.
        len_col : str
            length column name. Default is None.
        tiger_edges : bool
            using TIGER/Line edges file. Default is False.
        edge_subsets : list
            {type:{'col':column, 'oper':operator, 'val':value}}
            i.e. -- {'edge': {'col':'ROADFLG', 'val':'Y'}}
        mtfcc_split : str
            MTFCC codes for segments to weld and then split during the
            line splitting process. Default is None.
        mtfcc_intrst : str
            MTFCC codes for interstates. Default is None.
        mtfcc_ramp : str
            MTFCC codes for on ramps. Default is None.
        mtfcc_serv : str
            MTFCC codes for service drives. Default is None.
        mtfcc_split_grp : str
            after subseting this road type, group by this attribute
            before welding. Default is None.
        mtfcc_split_by : list
            MTFCC codes to eventually split the segments of
            `mtfcc_no_split` with. Default is None.
        skip_restr : bool
            skip re-welding restricted segments. Used when woring with
            TIGER/Lines. Default is False.
        tiger_roads : bool
            using TIGER/Line roads file. Default is False.
        calc_len : bool
            calculated length and add column. Default is False.
        record_components : bool
            record connected components in graph. This is used for
            teasing out the largests connected component.
            Default is False.
        largest_component : bool
            keep only the largest connected component in the graph.
            Default is False.
        record_geom : bool
            create associated between IDs and shapely geometries.
            Default is False.
        calc_stats : bool
            calculate network stats. Default is False.
        gen_matrix : bool
            calculate a cost matrix. Default is False.
        mtx_to_csv : str
            file path to save the cost matrix. Default is None.
        gen_paths : bool
            calculate shortest path trees. Default is False.
        paths_to_csv : str
            file path to save the shortest path trees. Default is None.
        gen_adjmtx : bool
            calculate node adjacency matrix. Default is False.
        adjmtx_to_csv : str
            file path to save the adjacency matrix. Default is None.
        algo : str
            shortest path algorithm. Default is None.
        def_graph_elems : bool
            define graph elements. Default is False.
        simplify : bool
            remove all non-articulation points from the network object.
            Default is False.
        save_full : bool
            save out the full network objects. Default is False.
        full_net_segms : str
            path and file name to save out the full network segments.
            Default is None.
        full_net_nodes : str
            path and file name to save out the full network nodes.
            Default is None.
        save_simplified : bool
            save out the simplified network objects. Default is False.
        simp_net_segms : str
            path and file name to save out the simplified network
            segments. Default is None.
        simp_net_nodes : str
            path and file name to save out the simplified network nodes.
            Default is None.
        remove_gdfs : bool
            remove dataframes from network object following  network
            simplification. Default is False.
        
        Methods : Attributes
        --------------------
        __init__ : segmdata, census_data
        build_network : --
        build_base : s_data, n_data, segm2xyid, node2xyid
        build_topology : segm2node, node2segment, segm2segm, node2node
        build_components : segm_cc, cc_lens, node_cc, longest_segm_cc,
            largest_segm_cc, largest_node_cc, n_edge_cc
        build_associations : s_ids, n_ids, n_segm, n_node, segm2len,
            network_length, node2degree, segm2tlid
        define_graph_elements : segm2elem, node2elem
        simplify_network : --
        granulate_network : ########??? add nodes for thiessen
        add_node : --
        add_edge : --
        adjacency_matrix : n2n_adjmtx
        network_cost_matrix : diameter, radius, d_net, d_euc, circuity,
            n2n_euclidean, n2n_algo, n2n_matrix, n2n_paths
        calc_net_stats : max_sinuosity, min_sinuosity,
            net_mean_sinuosity, net_std_sinuosity, max_node_degree,
            min_node_degree, mean_node_degree, std_node_degree, alpha,
            beta, gamma, eta, entropies_mtfcc, entropy_mtfcc,
            actual_object_sizes, actual_total_size
        sauce.setup_raw : raw_data_info
        sauce.ring_correction : corrected_rings
        sauce.line_splitter : lines_split
        sauce.seg_welder : welded_mls
        sauce.cleanse_supercycle : cleanse_cycles, scrubbed
        sauce.geom_assoc : segm2geom, node2geom
        sauce.coords_assoc : segm2coords, node2coords
        sauce.get_stats_frame : network_stats
        
        Examples
        --------
        >>> from spaghetti import spaghetti as spgh
        >>> net = spgh.SpaghettiNetwork()
        >>> print(net.network_stats)
        """
        
        
        if network_instance:
            self = network_instance
        else:
            self.tnid, self.tnidf, self.tnidt = tnid, tnidf, tnidt
            self.sid_name, self.nid_name = sid_name, nid_name
            self.len_col, self.geo_col = len_col, geo_col
            self.proj_init, self.proj_trans = proj_init, proj_trans
            self.proj_units = proj_units
            self.xyid = xyid
            self.inter = inter
            self.file_type = file_type
            self.mtfcc_types = mtfcc_types
            self.mtfcc_discard = mtfcc_discard
            self.tiger_edges = tiger_edges
            self.tiger_roads = tiger_roads
            self.discard_segs = discard_segs
            if self.tiger_edges or self.tiger_roads:
                self.census_data = True
            else:
                self.census_data = False
            
            if self.census_data:
                # TIGER variable attributes
                self.attr1 = attr1
                self.attr1rank = attr1rank
                self.attr2 = attr2
                self.study_area = study_area
                self.county = county
                self.state = state
                self.year = year
                self.place_time = place_time
                self.segmdata = segmdata
                if self.tiger_edges:
                    self.tlid = self.attr2
            
            # This reads in and prepares/cleans a segments geodataframe
            if not hasattr(segmdata, self.geo_col) and self.census_data:
                if self.tiger_edges:
                    self.edge_subsets = edge_subsets
                    self.mtfcc_split = mtfcc_split
                    self.mtfcc_intrst= mtfcc_intrst
                    self.mtfcc_ramp = mtfcc_ramp
                    self.mtfcc_serv = mtfcc_serv
                    self.mtfcc_split_grp = mtfcc_split_grp
                    self.mtfcc_split_by = mtfcc_split_by
                    self.skip_restr = skip_restr
                    
                    # fetch path to raw tiger edge data is available
                    # of local machine, otherwise download from
                    # https://www2.census.gov/
                    raw_file = sauce.get_raw_tiger_edges(self)
                elif self.tiger_roads:
                    pass
                
                # freshly cleaned segments geodataframe
                segmdata = sauce.tiger_netprep(self, in_file=raw_file,
                                               calc_len=calc_len)
            
            # build a network object from segments
            self.build_network(segmdata, record_components=record_components,
                               largest_component=largest_component,
                               record_geom=record_geom)
            
            if save_full:
                self.s_data.to_file(full_net_segms+self.file_type)
                self.n_data.to_file(full_net_nodes+self.file_type)
        
        # simplify the network
        if simplify:
            # create simplified segments geodataframe
            simplified_segms = self.simplify_network(self)
            # build a network object from simplified segments
            self.build_network(simplified_segms, record_geom=record_geom,
                               record_components=record_components,
                               largest_component=largest_component,
                               def_graph_elems=def_graph_elems)
            
            if save_simplified:
                self.s_data.to_file(simp_net_segms+self.file_type)
                self.n_data.to_file(simp_net_nodes+self.file_type)
        
        # create node to node adjacency matrix
        if gen_adjmtx:
            self.adjacency_matrix(adjmtx_to_csv)
        
        # Create and save out an all to all network node cost matrix
        if gen_matrix:
            
            if not algo:
                raise Exception('Set algorithm for cost matrix calculation.')
            
            self.n2n_algo = algo
            self.network_cost_matrix(mtx_to_csv=mtx_to_csv, gpths=gen_paths,
                                     paths_to_csv=paths_to_csv,
                                     calc_stats=calc_stats)
        # calculate descriptive network stats
        if calc_stats:
            self.calc_net_stats()
        
        if remove_gdfs:
            self.remove_frames()
    
    ###########################################################################
    ########################    end __init__    ###############################
    ###########################################################################
    
    def build_network(self, sdata, record_components=False, record_geom=False,
                      largest_component=False, def_graph_elems=False):
        """top-level method for full network object creation from a
        geopandas.GeoDataFrame of lines.
        
        Parameters
        ----------
        sdata : geopandas.GeoDataFrame
            segments data.
        record_components : bool
            find rooted connected components in the network (True),
            or ignore (False). Default is False.
        largest_component : bool
            keep only the largest connected compnent of the network
            (True), or keep all components (False). Default is False.
        record_geom : bool
            create an id to geometry lookup (True), or ignore (False).
            Default is False.
        def_graph_elems : bool
            define each element of the graph as either a branch
            [connected to two or more other elements], or a leaf
            [connected to only one other element] (True), or ignore
            (False). Default is False.
        """
        
        self.build_base(sdata)
        self.build_topology()
        if record_components:
            self.build_components(largest_cc=largest_component)
        self.build_associations(record_geom=record_geom)
        if def_graph_elems:
            self.define_graph_elements()
    
    
    def build_base(self, sdata):
        """extract nodes from segment endpoints and relate segments and
        nodes to a location ID (xyid)
        
        Parameters
        ----------
        sdata : geopandas.GeoDataFrame
            segments data.
        """
        
        # Instantiate segments dataframe as part of SpaghettiNetwork class
        self.s_data = sdata
        self.s_data.reset_index(drop=True, inplace=True)
        self.s_data = sauce.add_ids(self.s_data, id_name=self.sid_name)
        
        # create segment xyid
        self.segm2xyid = utils.generate_xyid(df=self.s_data, geom_type='segm',
                                             geo_col=self.geo_col)
        self.s_data = utils.fill_frame(self.s_data, idx=self.sid_name,
                                       col=self.xyid, data=self.segm2xyid)
        
        # Instantiate nodes dataframe as part of NetworkClass
        self.n_data = sauce.extract_nodes(self)
        self.n_data.reset_index(drop=True, inplace=True)
        
        # create permanent node xyid
        self.node2xyid = utils.generate_xyid(df=self.n_data, geom_type='node',
                                             geo_col=self.geo_col)
        self.n_data = utils.fill_frame(self.n_data, idx=self.nid_name,
                                       col=self.xyid, data=self.node2xyid)
    
    
    def build_topology(self):
        """relate all graph elements.
        """
        
        # Associate segments with neighboring nodes
        self.segm2node = sauce.associate(primary=self.segm2xyid,
                                         secondary=self.node2xyid,
                                         assoc='segm2node')
        
        # Associate nodes with neighboring segments
        self.node2segm = sauce.associate(primary=self.node2xyid,
                                         secondary=self.segm2xyid,
                                         assoc='node2segm')
        
        # Associate segments with neighboring segments
        self.segm2segm = sauce.get_neighbors(self.segm2node, self.node2segm,
                                             astype=list)
        
        # Associate nodes with neighboring nodes
        self.node2node = sauce.get_neighbors(self.node2segm, self.segm2node,
                                             astype=list)
        
        # 1. Catch cases with more than 2 neighboring nodes
        #    for a segment and throw up an error.
        # 2. Catch rings and add start and end node
        self = sauce.assert_2_neighs(self)
        
        # fill dataframe with seg2seg 
        self.s_data = utils.fill_frame(self.s_data, idx=self.sid_name,
                                       col='s_neigh', data=self.segm2segm)
        
        # fill dataframe with seg2node 
        self.s_data = utils.fill_frame(self.s_data, idx=self.sid_name,
                                       col='n_neigh', data=self.segm2node)
        
        # fill dataframe with node2seg 
        self.n_data = utils.fill_frame(self.n_data, idx=self.nid_name,
                                       col='s_neigh', data=self.node2segm)
        
        # fill dataframe with node2node 
        self.n_data = utils.fill_frame(self.n_data, idx=self.nid_name,
                                       col='n_neigh', data=self.node2node)
    
    
    def build_components(self, largest_cc=False):
        """find the rooted connected components of the graph (either
        largest or longest).  *** Must choose either largest or longest.
        If both `largest_cc` and `longest_cc` are `True`, `largest_cc`
        will be selected by default. ***
        
        Parameters
        ----------
        largest_cc : bool
            keep only the largest connected component (the most
            edges/nodes) in the graph. Default is False.
        """
        
        ### Segms -- Connected Components
        # -- Count
        self.segm_cc = sauce.get_roots(self.segm2segm)
        self.s_data = utils.fill_frame(self.s_data, idx=self.sid_name,
                                       col='CC', data=self.segm_cc)
        
        # -- Length
        # fill connected component len column in dataframe and return dict
        self.cc_lens = sauce.get_cc_len(self, len_col=self.len_col)
        
        ### Node -- Connected Components
        self.node_cc = sauce.get_roots(self.node2node)
        self.n_data = utils.fill_frame(self.n_data, idx=self.nid_name,
                                       col='CC', data=self.node_cc)
        
        # Extract largest CCs
        # largest CC by count (nodes & segments) -- and return small keys
        self.largest_segm_cc,\
        segm_smallkeys = sauce.get_largest_cc(self.segm_cc, smallkeys=True)
        self.largest_node_cc,\
        node_smallkeys = sauce.get_largest_cc(self.node_cc, smallkeys=True)
        
        # Keep on the largest connected component
        sauce.update_adj(self, segm_smallkeys, node_smallkeys)
        
        # Count connected components in network
        self.n_edge_cc = len(self.segm2segm)
    
    
    def build_associations(self, record_geom=False):
        """associate graph elements with geometries, coordinates,
        segment lengths, node degrees, and other information.
        
        Parameters
        ----------
        record_geom : bool
            create an id to geometry lookup (True).
        """
        
        if record_geom:
            sauce.geom_assoc(self)
        sauce.coords_assoc(self)
        
        # id lists
        self.s_ids = list(self.s_data[self.sid_name])
        self.n_ids = list(self.n_data[self.nid_name])
        
        # permanent segment count & node count
        self.n_segm, self.n_node = len(self.s_ids), len(self.n_ids)
        
        # Associate segments with length
        self.segm2len = sauce.xwalk(self.s_data, c1=self.sid_name,
                                    c2=self.len_col)
        
        # total length
        self.network_length = sum([v for (k,v) in self.segm2len])
        
        # Calculate degree for n_ids
        #   -- incident segs +1; incident loops +2
        self.node2degree = sauce.calc_valency(self, col='n_neigh')
        self.n_data['degree'] = [n2d[1][0] for n2d in self.node2degree]
        try:
            if self.tiger_edges:
                self.segm2tlid = sauce.xwalk(self.s_data, c1=self.sid_name,
                                             c2=self.tlid)
        except KeyError:
            pass
    
    
    def define_graph_elements(self):
        """define all segments and nodes as either a leaf (incident with
        one other element) or a branch (incident with more than one
        other grpah elemet). ** Only functional when the graph is
        comprised of one connected component.
        """
        
        self.segm2elem = sauce.branch_or_leaf(self, geom_type='segm')
        self.s_data = utils.fill_frame(self.s_data, idx=self.sid_name,
                                       col='graph_elem', data=self.segm2elem)
        self.node2elem = sauce.branch_or_leaf(self, geom_type='node')
        self.n_data = utils.fill_frame(self.n_data, idx=self.nid_name,
                                       col='graph_elem', data=self.node2elem)
    
    
    def simplify_network(self, save_full=False):
        """remove all non-articulation points in the network
        
        Parameters
        ----------
        save_full : bool
            save out full segments and nodes dataframe before replacing
            with simplified segments and nodes dataframe
        
        Returns
        -------
        simp_segms : geopandas.GeoDataFrame
            simplified segments dataframe.
        """
        
        # Create simplified road segments
        # (remove non-articulation points)
        simp_segms = sauce.simplify(self)
        
        # Reset index and SegIDX to match
        simp_segms.reset_index(drop=True, inplace=True)
        simp_segms = sauce.add_ids(simp_segms, id_name=self.sid_name)
        simp_segms = sauce.label_rings(simp_segms, geo_col=self.geo_col)
        fsave = self.inter+'Simplified_RingRoadsCorrected'
        simp_segms = sauce.ring_correction(self, simp_segms, save_keep=fsave)
        
        # add xyid
        segm2xyid = utils.generate_xyid(df=simp_segms, geom_type='segm',
                                        geo_col=self.geo_col)
        simp_segms = utils.fill_frame(simp_segms, col=self.xyid,
                                      data=segm2xyid)
        
        return simp_segms
    
    
    def calc_net_stats(self):
        """calculate network analyis descriptive statistics
        """
        
        # Calculate absolute shortest path along segments
        # Euclidean distance from vertex1 to vertex2
        self.s_data = sauce.euc_calc(self, col='euclid')
        
        # Calculate sinuosity for segments
        # curvilinear length / Euclidean Distance
        self.s_data['sinuosity'] = self.s_data[self.len_col]\
                                       /self.s_data['euclid']
        
        # networkSinuosity
        # set loop sinuosity (inf) to max nonloop sinuosity in dataset
        sinuosity = self.s_data['sinuosity'].copy()
        maxSin = sinuosity[sinuosity != np.inf].max()
        sinuosity = sinuosity.replace(to_replace=np.inf, value=maxSin)
        self.max_sinuosity = maxSin
        self.min_sinuosity = sinuosity.min()
        self.network_mean_sinuosity = sinuosity.mean()
        self.network_std_sinuosity = sinuosity.std()
        
        # node degree stats
        self.max_node_degree = self.n_data['degree'].max()
        self.min_node_degree = self.n_data['degree'].min()
        self.mean_node_degree = self.n_data['degree'].mean()
        self.std_node_degree = self.n_data['degree'].std()
        
        # network connectivity stats
        self.alpha = sauce.connectivity(self, measure='alpha')
        self.beta = sauce.connectivity(self, measure='beta')
        self.gamma = sauce.connectivity(self, measure='gamma')
        self.eta = sauce.connectivity(self, measure='eta')
        self.entropies_mtfcc = sauce.entropy(self) #return dict
        entropy = [v for k,v in list(self.entropies_mtfcc.items())]
        self.entropy_mtfcc = sum(entropy)*-1.
        
        # report object size
        self.actual_object_sizes = utils.ObjectSize(self)
        self.actual_total_size = self.actual_object_sizes.size_tot_all
        
        # create dataframe of descriptive network stats
        if hasattr(self, 'n2n_matrix'):
            sauce.get_stats_frame(self)
    
    
    def network_cost_matrix(self, calc_stats=False,  mtx_to_csv=None,
                            paths_to_csv=None, validate_symmetry=False,
                            gpths=False):
        """node to node cost matrix calculation with options for
        generating shortests paths along tree
        
        Parameters
        ----------
        calc_stats : bool
            calculate diameter and radius of the network
        validate_symmetry : bool
            validate matrix symmetry. Default is False.
        mtx_to_csv : str
            file path to save the cost matrix. Default is None.
        gpths : bool
            generate shortest paths tree. Default is False.
        paths_to_csv : str
            file path to save the paths. Default is None.
        """
        
        mtx, paths = sauce.shortest_path(self, gp=gpths)
        if calc_stats:
            
            # network diameter -- longest shortest path
            self.diameter = sauce._get_dia_rad(mtx, 'max')
            
            # network radius -- shortest shortest path
            self.radius = sauce._get_dia_rad(mtx, 'min')
            
            # circuity
            self.d_net = mtx.sum() # all network shortest paths
            coords = [v[0] for (k,v) in self.node2coords]
            self.n2n_euclidean = distance_matrix(coords, coords)
            
            # all euclidean shortest paths
            self.d_euc = self.n2n_euclidean.sum()
            self.circuity = self.d_net / self.d_euc
        
        if validate_symmetry:
            # validate symmetry
            if mtx[0][0] == 0.:
                if not sauce._check_symmetric(mtx, tol=1e-8):
                    raise Exception('all to all matrix is not symmetric')
        
        if mtx_to_csv:
            mtx_to_csv = mtx_to_csv + self.n2n_algo + '_CostMatrix'
            np.savetxt(mtx_to_csv+'.csv', mtx, delimiter=',')
        
        self.n2n_matrix = mtx
        
        if gpths:
            paths_to_csv = paths_to_csv + self.n2n_algo + '_AllShortestPaths'
            sauce._write_paths_csv(paths, paths_to_csv)
        
        if gpths:
            self.n2n_paths = paths
    
    
    def adjacency_matrix(self, adjmtx_to_csv):
        """Create a 1/0 adjacency matrix between all n_ids
        
        Parameters
        ----------
        adjmtx_to_csv   str
                        path to save adjacency matrix
        """
        
        self.n2n_adjmtx = np.zeros((self.n_node, self.n_node)).astype(int)
        
        for (i,js) in self.node2node:
            for j in js:
                self.n2n_adjmtx[i,j] = 1
        
        if adjmtx_to_csv:
            adjmtx = adjmtx_to_csv+'AdjacencyMatrix'
            np.savetxt(adjmtx+'.csv', self.n2n_adjmtx, fmt='%i', delimiter=',')
    
    
    def remove_frames(self):
        """remove dataframes from network object
        """
        del self.s_data
        del self.n_data


class SpaghettiPointPattern():
    """
    keeping as a class for now to mimic pysal
    ...but is it really mimicing now?
    --> Network K , Cross
    --> Network Nearest Neighbor, Cross
    """
    
    
    def __init__(self, df=None, df_name=None, df_key=None, df_pop=None,
                 net_segms=None, net=None, simulated=False,
                 kd_tree=None, k=5, tol=.01, snap_to='segments',
                 no_pop=['FireStations', 'FireStationsSynthetic']):
        """
        Parameters
        ----------
        df : geopandas.GeoDataFrame
            observation points dataframe. Default is None.
        df_name : str
            dataframe name. Default is None.
        df_key : str/int
            dataframe key column name. Default is None.
        net : spaghetti.SpaghettiNetwork
            network object. Default is None.
        net_segms : geopandas.GeoDataFrame
            network segments dataframe. Default is None.
        simulated : bool
            empir. or sim. points along network segments.
            Default is False.
        kd_tree : scipy.spatial.kdtree.KDTree
            all network nodes lookup
        k : int
            number of nearest neighbors to query. Default is 5.
        tol : float
            snapping to line tolerance. Default is .01.
        snap_to : str
            snap points to either segments of nodes.
            Default is 'segments'.
        no_pop : list
            observations that do not include a population measure.
            Default is `['FireStations', 'FireStationsSynthetic']`.
        
        Methods : Attributes
        --------------------
        study_area : str
            study area within county
        sid_name : str
            segment id column name.
        proj_init : int
            initial projection
        xyid : str
            combined x-coord + y-coords string ID
        geo_col : str
            geometry column name
        obs2coords : list
            observation index and attribute id lookup of coordinates.
        snapped_points : geopandas.GeoDataFrame
            snapped point representation
        obs2segm : dict
            observation id (key) to segment id
        """
        
        self.study_area = net.study_area
        self.sid_name = net.sid_name
        self.proj_init = net.proj_trans
        self.df = df
        self.df_name = df_name
        self.df_key = df_key
        self.xyid = net.xyid
        self.geo_col = net.geo_col
        self.k = k
        self.kd_tree = kd_tree
        self.tol = tol
        self.snap_to = snap_to
        
        # create observation to coordinate xwalk
        self.obs2coords = sauce.get_obs2coords(self)
        # snap points and return dataframe
        self.snapped_points = sauce.snap_to_nearest(self, net=net,
                                                    sframe=net_segms)
        # set dataframe projection
        self.snapped_points = utils.set_crs(self.snapped_points,
                                            proj_init=self.proj_init)
        
        self.obs2segm = dict(zip(self.snapped_points[self.df_key],
                                 self.snapped_points['assoc_segm']))
        
        if not self.df_name in no_pop:
            
            try:
                self.snapped_points[df_pop] = self.df[df_pop]
            except KeyError:
                try:
                    df_pop = 'POP100_syn'
                    self.snapped_points[df_pop] = self.df[df_pop]
                except KeyError:
                    df_pop = 'synth_pop'
                    self.snapped_points[df_pop] = self.df[df_pop]
            
            #create a segment-to-population tracker
            # this will vary with each different method employed
            self.segm2pop = {seg: self.snapped_points.loc[\
                                  (self.snapped_points['assoc_segm'] == seg),
                                                       df_pop].sum()
                                                       for seg in net.s_ids}


def remove_restricted_segms(df, net, restr=None, col=None):
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
    df : geopandas.GeoDataFrame
        network segments of unrestricted edges dataframe
    net : spaghetti.SpaghettiNetwork
    """
    
    df, net = sauce.remove_restricted(df, net, restr=restr, col=col)
    
    return df, net


def build_net_nodes_kdtree(nodes, only_coords=False):
    """build a kdtree from the network node coords for
    point pattern lookup
    
    Parameters
    ----------
    node : either nodes2coords dictionary or coords list
    only_coords : bool
    
    Returns
    -------
    net_nodes_kdtree : scipy.spatial.kdtree.KDTree
        all network nodes lookup
    """
    
    if only_coords:
        return cg.KDTree(nodes)
    
    else:
        if type(nodes) is dict:
            node_coords = [coords[0] for node,coords in list(nodes.items())]
        elif type(nodes) is list:
            node_coords = [coords[0] for (node,coords) in nodes]
    
    net_nodes_kdtree = cg.KDTree(node_coords)
    
    return net_nodes_kdtree


def obs2obs_cost_matrix(orig, dest=None, dist_col=None, symmetric=False,
                        network_matrix=None, from_nodes=False, wsnap_dist=None,
                        assoc_col=None, distance_type=None, xyid=None,
                        numeric_cols=None):
    """calculate a cost matrix from (n) observations to (m) observations
    
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
    
    n2m_matrix = sauce.obs2obs_costs(orig, dest=dest, symmetric=symmetric,
                                     network_matrix=network_matrix, xyid=xyid,
                                     from_nodes=from_nodes,
                                     wsnap_dist=wsnap_dist,
                                     distance_type=distance_type,
                                     assoc_col=assoc_col,
                                     numeric_cols=numeric_cols)
    
    return n2m_matrix


def obs2obs_nearestneighbor(matrix, orig_ids=None, dest_ids=None,
                            symmetric=False, keep_zero_dist=True):
    """Compute the interpattern nearest neighbor distances or the
    intrapattern nearest neighbor distances between a source
    pattern and a destination pattern.
    
    Parameters
    ----------
    matrix : numpy.ndarray
        nXm cost matrix
    orig_ids : list
        'descriptive' observation ids (not simple index values)
    dest_ids : list
        'descriptive' observation ids (not simple index values)
    symmetric : bool
        origin to origin (True), origin to destination (False).
    keep_zero_dist : bool
        Include zero values in minimum distance (True) or exclude
        (False). Default is True. If the source pattern is the same as
        the destination pattern the diagonal is filled with nans
    
    Returns
    -------
    nearest : dict
        origin (descriptive id, matrix index) key and destination
        ((descriptive id list, matrix index list), distance)
        for the value
    """
    
    # (for source-to-source patterns) if zero-distance neighbors should
    # be ignored, convert keep the diagonal as 0.0 and take the minimum
    # distance neighbor(s) that is/are not 0.0 distance.
    if keep_zero_dist and symmetric:
        # (for source-to-source patterns) if zero-distance neighbors are
        # desired, set the diagonal to NaN and take the minimum distance
        # neighbor(s), which may include zero distance neighors.
        np.fill_diagonal(matrix, np.nan)
    
    nearest = {}
    for src_idx, orig_id in enumerate(orig_ids):
        if keep_zero_dist:
            val = np.nanmin(matrix[src_idx,:])
        else:
            val = np.min(matrix[src_idx,:][np.nonzero(matrix[src_idx,:])])
        
        # nearest destination (may be more than one if equal distance)
        dst_idxs = np.where(matrix[src_idx,:] == val)[0].tolist()
        
        # fetch descriptive id
        dest_desc_ids = [dest_ids[idx] for idx in dst_idxs]
        
        # set nearest with origin (descriptive id, matrix index) key and
        # destination ((descriptive id list, matrix index list),
        # distance) for the value
        nearest[(orig_id, src_idx)] = ((dest_desc_ids,dst_idxs), val)
        
    return nearest


def synthetic_grid(bounds=None, n_hori_lines=None, n_vert_lines=None,
                   withbox=False, as_polys=False):
    """generate a graph theoretic lattice. set kwargs `withbox`
    and `as_polys` to create grid of polygons.
     
    Parameters
    ----------
    bounds : list
        area bounds in the form of [x1,y1,x2,y2].
    n_hori_lines : int
        Count of horizontal lines. Default is None.
    n_vert_lines : int
        Count of vertical lines. Default is None.
    withbox : bool
        Include outer bounding box segments. Default is False.
    as_polys : bool
        polygonize network blocks. Default is False.
   
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        either LineStrings or Polygons
    """
    
    # bounding box line lengths
    h_length, v_length = bounds[2] - bounds[0],\
                         bounds[3] - bounds[1]

    # horizontal and vertical increments
    h_incr, v_incr = h_length/float(n_hori_lines+1),\
                     v_length/float(n_vert_lines+1)

    # define the horizontal and vertical space
    hspace = [h_incr*slot for slot in range(n_vert_lines+2)]
    vspace = [v_incr*slot for slot in range(n_hori_lines+2)]

    # get vertical and horizontal lines
    horis = sauce._get_lines(hspace, vspace, withbox, bounds)
    verts = sauce._get_lines(hspace, vspace, withbox, bounds, hori=False)

    # combine into one list and 
    geoms = verts + horis
    
    if as_polys and withbox:
        geoms = sauce._polygonize(geoms)
    
    # instantiate dataframe
    gdf = gpd.GeoDataFrame(geometry=geoms)
    
    return gdf


def dump_pickled(obj, pickle_jar, pickle_name='Network'):
    """pickle an object.
       if on MacOSX work in chunks due to bug for reading/writing over
       2GB. https://bugs.python.org/issue24658
    
    Parameters
    ----------
    obj : any object
    pickle_jar : str
        path to the pickled object
    pickle_name : str
        name of pickled object. Default is 'Network'.
    """
   
    out_file = pickle_jar + pickle_name + '.pkl'
    sys_os = platform.system()
    bytes_maximum = 2**31 - 1
    obj_in_bytes_out = pickle.dumps(obj, protocol=4)
    
    if (len(obj_in_bytes_out) <= bytes_maximum and sys_os == 'Darwin')\
    or sys_os != 'Darwin':
        pickled_obj = open(out_file, 'wb')
        pickle.dump(obj, pickled_obj, protocol=pickle.HIGHEST_PROTOCOL)
        pickled_obj.close()
    
    else:
        with open(out_file, 'wb') as obj_out:
            for idx in range(0, len(obj_in_bytes_out), bytes_maximum):
                obj_out.write(obj_in_bytes_out[idx:idx+bytes_maximum])


def load_pickled(pickle_jar, pickle_name='Network'):
    """load a pickled object
        if on MacOSX work in chunks due to bug for reading/writing over
        2GB. https://bugs.python.org/issue24658
    
    Parameters
    ----------
    pickle_jar : str
        path to the pickled object
    pickle_name : str
        name of pickled object. Default is 'Network'.
    
    Returns
    -------
    unpickled spaghetti.Spaghetti network object
    """
    
    in_file = pickle_jar + pickle_name + '.pkl'
    if not os.path.exists(in_file):
        raise FileNotFoundError
    
    sys_os = platform.system()
    bytes_maximum = 2**31 - 1
    obj_in_bytes_in = bytearray(0)
    obj_size = os.path.getsize(in_file)
    if (obj_size <= bytes_maximum and sys_os == 'Darwin')\
    or sys_os != 'Darwin':
        pickled_obj = open(in_file, 'rb')
        return pickle.load(pickled_obj)
    
    else:
        with open(in_file, 'rb') as file_in:
            for _ in range(0, obj_size, bytes_maximum):
                obj_in_bytes_in += file_in.read(bytes_maximum)
        return pickle.loads(obj_in_bytes_in)

