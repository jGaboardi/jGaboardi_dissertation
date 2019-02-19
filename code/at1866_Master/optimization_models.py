"""Optimal facility Location Modeling
"""

# standard library imports
from collections import OrderedDict
import time

# non-standard library imports
import geopandas as gpd
from libpysal import cg
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
from PIL import Image
import pulp
import seaborn as sns
from shapely.geometry import Point, LineString, MultiLineString


class FacilityLocationModel:
    """Solve a facility locatiom optimization model
    
    Parameters 
    ----------
    locmod : str
        location model name; must also be defined as a class method.
    name : str
        full model name including spatial extent, clients type,
        facility type, parameters, etc. Default is None.
    cij : numpy.ndarray
        cost matrix from origins (index of i) to destination
        (index of j). Default is None.
    ai : numpy.ndarray
        Client weight vector. Default is None.
    origin_service : bool
        Origins are service facilities and destinations are
        clients [True]. origins are clients and destinations are
        service facilities [False]. Default is True.
    s : float
        service radius. Default is None.
    p : int
        Density of facilities to site. Default is None.
    write_lp : str
        file name (and path) of the LP file to write out.
    solver : str
        colloquial name of solver to be used
    print_sol : bool
        print select results. Default is True.
    cli_dvs : iterable
        client decision variable names. Default is None.
    fac_dvs : iterable
        service decision variable names. Default is None.
    frac_gap : float
        CBC param. fractional gap tolerence for solution.
    time_limit : {int, float}
        CBC param. time limit for solution in minutes, which
        is then convertered to seconds when passed to CBC.
    
    Methods
    -------
    build_lscp : build location set covering problem
    build_pmp : build p-median problem
    build_pcp : build p-center problem
    build_mclp : build maximal covering location problem
    problem_instance : create a problem instance
    add_vars : add variables to a model
    add_constrs : add contraints to a model
    add_obj : add an objective function to a model
    optimize : solve a model
    record_decisions : record optimal decision variables
    non_obj_vals : record non-objective values stats
        (eg. percent covered)
    print_results : print selected results
    create_cli2fac_sp : add shortest path lookup between orgin
        and destination point patterns.
    
    Attributes
    ----------
    problem : pulp.pulp.LpProblem
        problem/model instance
    n_cli : int
        total client sites
    r_cli : range
        iterable of client sites
    n_fac : int
        total candidate facility sites
    r_fac : range
        iterable of candidate facility sites
    aij : numpy.ndarray
        binary coverage matrix from cij (within s service radius)
    sij : numpy.ndarray
        demand weighted cost matrix as (ai * cij).
    cli_decisions : list, np.array
        client-to-service decision variable matrix used in pmp and pcp.
    fac_vars : dict
        facility decision variables
    cli_vars : dict
        client decision variables
    W : pulp.pulp.LpVariable
        minimized maximum variable in the p-center problem formulation
    build_time : float
        minutes to build the model
    solve_time : float
        minutes to solve the model
    record_decisions_time : 
        minutes to record model decision variables
    total_time : float
        total run minutes
    lp_formulation : str
        linear programming formulation of the model
    solve_minutes : float
        solve time in minutes
    obj_val : int or float
        model objective value
    fac2cli : dict
        facility to client relationship lookup
    cli2fac : dict
        client to facility relationship lookup
    fac2iloc : dict
        facility to dataframe index location lookup
    n_cli_uncov : int
        count of client location outside the service radius  
    cli2ncov : dict
        client to covered by count lookup
    ncov2ncli : dict
        covered by count to client count lookup
    mean_dist :
        mean distance per person to the assigned facility
    perc_served : 
        percentage of weighted clients covered in `s`
    cli2fac_sp : dict
        snapped point pattern OD network node lookup
    gap : {int, float, str}
        when a problem proven optimal this is an integer (0) or
        a very small float. Else it is a string.
    solver : class
        instance of solver
    solver_name : str
        internal pulp name of solver
    """
    
    def __init__(self, locmod, name=None, ai=None, cij=None, s=None,
                 p=None, write_lp=None, print_sol=True, cli_dvs=None,
                 fac_dvs=None, frac_gap=None, time_limit=None, solver='CBC',
                 origin_service=True):
        
        # Set model information
        self.locmod = locmod
        if name:
            self.name = name
        else:
            self.name = self.locmod
        
        # set CBC params
        self.frac_gap = frac_gap
        self.time_limit = time_limit * 60. if time_limit else None
        
        # Set parameters and indices
        # facility parameter
        self.p = p
        
        # client and facility count
        self.cij = cij
        if origin_service:
            self.n_cli, self.n_fac = cij.shape[1], cij.shape[0]
        else:
            self.n_cli, self.n_fac = cij.shape[0], cij.shape[1]
        
        # client and facility range
        self.r_cli = range(self.n_cli)
        self.r_fac = range(self.n_fac)
        
        # demand parameter
        if ai is not None:
            self.ai = ai.values
            self.ai_sum = ai.sum()
            
            # weighted demand
            self.sij = np.multiply(self.ai, self.cij)
        
        # if the model has a service radius parameter
        self.s = s
        if self.s:
            # binary coverage matrix from cij
            self.aij = np.zeros(self.cij.shape)
            self.aij[self.cij <= self.s] = 1.
        
        # if client and service decision variable names are passed
        # in the class, set them as attributes here.s
        if fac_dvs is not None:
            self.fac_dvs = list(fac_dvs)
        if cli_dvs is not None:
            self.cli_dvs = list(cli_dvs)
        
        # Set decision variables, constraints, and objective function
        try:
            self.build_time = getattr(self, 'build_'+self.locmod)()
        except AttributeError:
            raise AttributeError(self.locmod,
                                 'not a defined location model.')
        
        # solve
        self.solve_time = self.optimize(solver=solver, write_lp=write_lp)
        
        # records seleted decision variables
        self.record_decisions_time = self.record_decisions()
        
        # record non-objective values stats (eg. percent covered)
        if self.locmod in ['mclp', 'pmp']:
            self.non_obj_vals()
        self.total_time = self.build_time\
                          + self.solve_time\
                          + self.record_decisions_time
        
        # print results
        if print_sol:
            self.print_results()
    
    
    def timer(timed_function):
        """ decoration wrapper to time functions
        """
        
        def timer(*args, **kwargs):
            
            t1 = time.time()
            timed_function(*args, **kwargs)
            t2 = time.time() - t1
            timer = t2/60.
            
            return timer
        
        return timer
    
    
    @timer
    def build_lscp(self):
        """ Integer programming formulation of the Location Set
        Covering Problem.
        
        Originally Published:
            Toregas, C. and ReVelle, Charles. 1972.
            Optimal Location Under Time or Distance Constraints.
            Papers of the Regional Science Association.
            28(1):133 - 144.
        """
        
        # Problem Instance
        self.problem_instance()
        
        # Decision Variables
        self.add_vars()
        
        # Constraints
        self.add_constrs(constr=1) # set coverage constraints
        
        # Objective Function
        self.add_obj()
    
    
    @timer
    def build_mclp(self):
        """Integer programming formulation of the Maximal
        Covering Location Problem.
        
        Originally Published:
            Church, R. L and C. ReVelle. 1974. The Maximal
            Covering Location Problem. Papers of the Regional
            Science Association. 32:101-18.
        """
        
        # Problem Instance
        self.problem_instance()
        
        # Decision Variables
        self.add_vars()
        
        # Constraints
        self.add_constrs(constr=3) # facility constraint
        self.add_constrs(constr=6) # maximal coverage constraints
        
        # Objective Function
        self.add_obj()
    
    
    @timer
    def build_pmp(self):
        """Integer programming formulation of the p-median Problem.
        
        Originally Published:
            S. L. Hakimi. 1964. Optimum Locations of Switching
            Centers and the Absolute Centers and Medians of a Graph.
            Operations Research. 12 (3):450-459.
        
        Adapted from:
                -1-
            ReVelle, C.S. and Swain, R.W. 1970. Central facilities
            location. Geographical Analysis. 2(1), 30-42.
                -2-
            Toregas, C., Swain, R., ReVelle, C., Bergman, L. 1971.
            The Location of Emergency Service Facilities. Operations
            Research. 19 (6), 1363-1373.
                - 3 -
            Daskin, M. (1995). Network and discrete location:
            Models, algorithms, and applications. New York:
            John Wiley and Sons, Inc.
        """
        
        # Problem Instance
        self.problem_instance()
        
        # Decision Variables
        self.add_vars()
        
        # Constraints
        self.add_constrs(constr=2) # assignment constraints
        self.add_constrs(constr=3) # facility constraint
        self.add_constrs(constr=4) # opening constraints
        
        # Objective Function
        self.add_obj()
    
    
    @timer
    def build_pcp(self):
        """Integer programming formulation of the p-center Problem.
        
        Originally Published:
            S. L. Hakimi. 1964. Optimum Locations of Switching
            Centers and the Absolute Centers and Medians of a Graph.
            Operations Research. 12 (3):450-459.
        
        Adapted from:
            Daskin, M. (1995). Network and discrete location:
            Models, algorithms, and applications. New York:
            John Wiley and Sons, Inc.
        """
        
        # Problem Instance
        self.problem_instance()
        
        # Decision Variables
        self.add_vars()
        
        # Constraints
        self.add_constrs(constr=2) # assignment constraints
        self.add_constrs(constr=3) # facility constraint
        self.add_constrs(constr=4) # opening constraints
        self.add_constrs(constr=5) # minimized maximum constraints
        
        # Objective Function
        self.add_obj()
    
    
    def problem_instance(self):
        """Create a problem instance and set objective sense.
        """
        
        # set sense
        if self.locmod == 'mclp':
            sense = pulp.LpMaximize
        else:
            sense = pulp.LpMinimize
        
        # create problem instance
        self.problem = pulp.LpProblem(self.name, sense)
    
    
    def add_vars(self):
        """Add variables to a model.
        """
        
        # facility decision variables
        # variable names
        if not hasattr(self, 'fac_dvs'):
            self.fac_dvs = ['y%s' % i for i in self.r_fac]
        
        # pulp.LpVariable objects
        self.fac_vars = [pulp.LpVariable(var_name, cat='Binary')\
                                            for var_name in self.fac_dvs]
        
        # add variables
        self.problem.addVariables(v for v in self.fac_vars)
        
        # client decision variables
        # variable names
        if not hasattr(self, 'cli_dvs'):
            self.cli_dvs = ['x%s' % j for j in self.r_cli]
        
        if self.locmod == 'mclp':
            # pulp.LpVariable objects
            self.cli_vars = [pulp.LpVariable(var_name, cat='Binary')\
                                    for var_name in self.cli_dvs]
            # add variables
            self.problem.addVariables(v for v in self.cli_vars)
        
        if self.locmod in ['pmp', 'pcp']:
            self.cli_decisions = [['x%s_%s' % (i,j) for i in self.r_fac]\
                                                    for j in self.r_cli]\
            # pulp.LpVariable objects
            self.cli_vars = [[pulp.LpVariable(var_name, cat='Binary')\
                                             for var_name in var_row]\
                                             for var_row in self.cli_decisions]
            
            # add variables
            self.problem.addVariables(var for var_row in self.cli_vars\
                                          for var in var_row)
        
        # minimized maximum variable
        if self.locmod == 'pcp':
            self.W = pulp.LpVariable('W', cat='Continuous')
            self.problem.addVariable(self.W)
        
        self.fac_vars = np.array(self.fac_vars)
        if hasattr(self, 'cli_vars'):
            self.cli_vars = np.array(self.cli_vars)
    
    
    def add_constrs(self, constr=None):
        """ Add constraints to a model.
        (1) set coverage constraints
                y1 + x2 >= 1
                x1 + x3 >= 1
                x2 >= 1
        (2) assignment constraints
                x1_1 + x1_2 + x1_3 = 1
        (3) facility constraints
                y1 + y2 + y3 = p
        (4) opening constraints
                - x1_1 + y1 >= 0
                - x2_1 + y1 >= 0
                - x3_1 + y1 >= 0
        (5) minimax constraints
                cost1_1*x1_1 + cost1_2*x1_2 + cost1_3*x1_3 - W <= 0
        (6) maximal coverage constraints
                - x1 + y1 + y3 >= 0
                - x2 + y4 >= 0
        
        Parameters
        ----------
        constr : int {1, 2, 3, 4, 5, 6}
            Contraint type to add to model. See above for explanation.
            Default is None.
        """
        
        # 1 - set covering constraints
        if constr == 1:
            constrs = [pulp.lpSum(self.aij[i,j] * self.fac_vars[i]\
                                       for i in self.r_fac) >= 1\
                                       for j in self.r_cli]
        
        # 2 - assignment constraints
        elif constr == 2:
            constrs = [pulp.lpSum(self.cli_vars[j][i]
                                  for i in self.r_fac) == 1\
                                  for j in self.r_cli]
        
        # 3 - facility constraint
        elif constr == 3:
            constrs = [pulp.lpSum(self.fac_vars[i]\
                                               for i in self.r_fac) == self.p]
        
        # 4 - opening constraints
        elif constr == 4:
            constrs = [self.fac_vars[i] >= self.cli_vars[j][i]\
                       for j in self.r_cli for i in self.r_fac]
        
        # 5 - minimax constraints
        elif constr == 5:
            constrs = [pulp.lpSum(self.cij[i,j]\
                                  * self.cli_vars.T[i][j]
                                  for i in self.r_fac) <= self.W\
                                  for j in self.r_cli]
        
        # 6 - max coverage constraints
        elif constr == 6:
            constrs = [pulp.lpSum(self.aij[i,j]\
                                  * self.fac_vars[i]\
                                  for i in self.r_fac)\
                                  >= self.cli_vars[j]
                                  for j in self.r_cli]
        
        [self.problem.addConstraint(c) for c in constrs]
    
    
    def add_obj(self):
        """ Add an objective function to a model.
        """
        
        if self.locmod == 'lscp':
            self.problem.setObjective(pulp.lpSum(self.fac_vars[j]\
                                                 for j in self.r_fac))
        
        elif self.locmod == 'mclp':
            flat_ai = self.ai.flatten()
            self.problem.setObjective(pulp.lpSum(flat_ai[j] * self.cli_vars[j]
                                                 for j in self.r_cli))
        
        elif self.locmod == 'pmp':
            self.problem.setObjective(pulp.lpSum(self.sij[i,j]\
                                      * self.cli_vars.T[i][j]\
                                      for j in self.r_cli\
                                      for i in self.r_fac))
        
        elif self.locmod == 'pcp':
            self.problem.setObjective(self.W)
    
    
    @timer
    def optimize(self, solver=None, write_lp=None):
        """ Solve the model.
        
        Parameters
        ----------
        write_lp : str
            write out the linear programming formulation
        """
        cbc_variants = ['PULP', 'CBC', 'COIN', 'COIN-CBC', 'PULP_CBC_CMD']
        
        if solver.upper() in cbc_variants:
            self.solver = pulp.PULP_CBC_CMD
            self.solver_name = solver
        
        self.problem.solve(self.solver(maxSeconds=self.time_limit,
                                             fracGap=self.frac_gap))
        
        # linear programming formulation
        if write_lp:
            self.problem.writeLP(write_lp+'.lp')
        
        self.obj_val = self.problem.objective.value()
    
    
    @timer
    def record_decisions(self):
        """record decision variable relationship
        folowing optimization.
        """
        
        # facility-to-dataframe index location lookup
        self.fac2iloc = {v.name: idx for idx, v\
                                     in enumerate(self.fac_vars)}
        
        # client-to-dataframe index location lookup
        self.cli2iloc = {}
        
        # facility-to-client lookup
        self._fac2cli()
        
        # client-to-facility lookup
        self._cli2fac()
        
        # count of uncovered clients
        self.n_cli_uncov = self.n_cli - len(self.cli2iloc.keys())
        
        # count of clients covered by n facilities
        if self.locmod in ['lscp', 'mclp']:
            self._cli2cov()
    
    
    def _fac2cli(self):
        """internal facility-to-client lookup function
        """
        
        self.fac2cli = {}
        for i in self.r_fac:
            if self.fac_vars[i].value() > 0:
                jvar = str(self.fac_dvs[i])
                self.fac2cli[jvar] = []
                for j in self.r_cli:
                    ivar = None
                    if self.locmod == 'lscp':
                        if self.aij[i,j] > 0:
                            ivar = str(self.cli_dvs[j])
                            self.fac2cli[jvar].append(ivar)
                    elif self.locmod == 'mclp':
                        if self.cli_vars[j].value() > 0:
                            if self.aij[i,j] > 0:
                                ivar = str(self.cli_dvs[j])
                                self.fac2cli[jvar].append(ivar)
                    else:
                        if self.cli_vars[j][i].value() > 0:
                            ivar = str(self.cli_dvs[j])
                            self.fac2cli[jvar].append(ivar)
                    if ivar:
                        self.cli2iloc[ivar] = j
    
    
    def _cli2fac(self):
        """internal client-to-facility lookup function
        """
        
        self.cli2fac = {}
        
        for cv in list(self.cli2iloc.keys()):
            self.cli2fac[cv] = []
            for k,v in self.fac2cli.items():
                if cv in v:
                    self.cli2fac[cv].append(k)
        
        
    def _cli2cov(self):
        """internal client coverage function
        """
        
        # client decision variable: coveraged by count
        self.cli2ncov = {}
        
        # tease out how many facilities cover each client
        for c, fs in self.cli2fac.items():
            self.cli2ncov[c] = len(fs)
        
        # determine coverage count of each facility
        try:
            most_coverage = max(self.cli2ncov.values())
        except ValueError:
            raise Exception('Service radius of %s covers no clients' % self.s)
        
        # covered client total by facility total
        # "3 clients are covered by 2 facilities"
        self.ncov2ncli = {}
        
        for cov_count in range(most_coverage+1):
            
            # uncovered clients
            if cov_count == 0:
                self.ncov2ncli[cov_count] = self.n_cli_uncov
                continue
            
            if not cov_count in list(self.cli2ncov.keys()):
                self.ncov2ncli[cov_count] = 0
            
            for c, ncov in self.cli2ncov.items():
                if ncov >= cov_count:
                    self.ncov2ncli[cov_count] += 1
    
    
    def non_obj_vals(self):
        """non-objective value metrics
        """
        
        mean_stat = self.obj_val/float(self.ai_sum)
        
        if self.locmod == "pmp":
            self.mean_dist = mean_stat
        if self.locmod == "mclp":
            self.perc_served = mean_stat * 100.
    
    
    def print_results(self):
        """print select results
        """
        
        print('\nModel:\t%s' % self.name)
        
        # issue statement if model not optimal
        self.status_int = self.problem.status
        self.status = pulp.LpStatus[self.status_int]
        if self.status_int != 1:
            warning = '  !!! WARNING!!!  '
            print('see [' + warning + '] below')
            print('\t%s model is %s.' % (self.name, self.status))
            params = []
            for param in ['p','s']:
                val = getattr(self, param)
                if val:
                    params.append(param+'=%s'%val)
            
            print('\t\tparameters:', ', '.join(params))
            print('\tThe objective function has been relaxed.\n')
        
        # coverage values
        if self.locmod in ['lscp', 'mclp']:
            for ncov, ncli in self.ncov2ncli.items():
                if ncov == 0:
                    cov_msg = 'client locations are not covered'
                    if self.locmod == 'lscp' and self.n_cli_uncov > 0:
                        warning = ' <------------' + warning
                        cov_msg = cov_msg.upper() + warning
                    print('--- %i' % self.n_cli_uncov, cov_msg)
                else:
                    if ncov == 1:
                        sp = 'y'
                    else:
                        sp = 'ies'
                    print('--- %i client locations are covered' % ncli,
                          'by %i' % ncov, 'facilit'+sp)
        
        # solve time and objective value
        if self.locmod == 'lscp':
            if self.status_int == 1:
                cov_type = ' total '
            else:
                cov_type = ' partial '
            u1 = 'facilities needed for' + cov_type + 'coverage within a '
            u2 = '%s meter service radius' % self.s
        else:
            if self.locmod == 'pmp':
                u1 = 'total weighted distance with '
            if self.locmod == 'pcp':
                u1 = 'worst case distance with '
            if self.locmod == 'mclp':
                u1 = 'residents within %s meters of '% self.s
            u2 = '%i selected facilities' % self.p
        units = u1 + u2
        
        print('Obj. Value:\t%s %s' % (round(self.obj_val, 3), units))
        rounded_time = round(self.solve_time, 3)
        print('Build Time:\t%s minutes' % round(self.build_time, 3))
        print('Solve Time:\t%s minutes' % rounded_time)
        print('Total Time:\t%s minutes' % round(self.total_time, 3))
        if self.time_limit:
            print('Time Limit:\t%s minutes' % str(self.time_limit/60.))
        try:
            self.gap = self.problem.infeasibilityGap()
        except ValueError:
            self.gap = 'sub-optimal solution'
        print('Infeas. Gap:\t%s' % self.gap)
        if self.locmod == 'pmp':
            print('Mean distance per',
                  'person: %f' % self.mean_dist)
        if self.locmod == 'mclp':
            print('Percent of %i' % self.ai_sum,
                  '(population) covered: %f' % self.perc_served)
    
    
    def create_cli2fac_sp(self, tree):
        """add shortest path lookup between orig. and dest. point patterns.
        
        Parameters
        ----------
        tree : dict
            origin node to desitination network node
        index look up keyed on OD array iloc.
        """
        
        self.cli2fac_sp = {}
        for fac, clis in self.fac2cli.items():
            for cli in clis:
                self.cli2fac_sp[fac,cli] = []
                facloc = self.fac2iloc[fac]
                cliloc = self.cli2iloc[cli]
                try:
                    self.cli2fac_sp[fac,cli] = tree[(cliloc, facloc)]
                except:
                    self.cli2fac_sp[fac,cli] = ('snap', 'snap')


def element_calculator(locmod, rows, cols=None, print_count=False):
    """print the number of constraints and variables in the program.
    
    Parameters
    ----------
    locmod : str
        model name
    rows : int
        client locations
    cols : int
        service locations
    
    Returns
    -------
    n_constraints : int
        number of constraints
    n_variables : int
        number of constraints
    """
   
    n_constraints = 0
    n_variables = 0
    
    if locmod == 'lscp':
        n_constraints += cols
        n_variables += rows
    
    if locmod in ['mclp', 'pmp', 'pcp']:
        facility = 1
        n_constraints += facility
        n_variables += rows
        
    if locmod in ['pmp', 'pcp']:
        assignment = cols
        opening = cols * rows
        n_constraints += assignment + opening
        n_variables += cols * rows
    
    if locmod == 'pcp':
        minimax = cols
        w = 1
        n_constraints += minimax
        n_variables += w
        
    if locmod == 'mclp':
        max_coverage = cols
        n_constraints += max_coverage
        n_variables += max_coverage
        
    if print_count:
        print('n Constraints:', n_constraints, 'element_calculator')
        print('n Variables:', n_variables, 'element_calculator')
        print('\n')
    
    return n_constraints, n_variables


def add_results(model, cli_df, fac_df,
                print_solution=False, sol_col='solution'):
    """Add decision variable relationships to a dataframe.
    
    Parameters
    ----------
    model : at1866_Master.FacilityLocationModel
    cli_df : geopandas.GeoDataFrame
        GeoDataFrame of client locations
    fac_df : geopandas.GeoDataFrame
        GeoDataFrame of facility locations
    print_solution : bool
        print out solution decision variables. Default is False.
    sol_col : str
        solution column name. Default is 'solution'.
    
    Returns
    -------
    cli_df : geopandas.GeoDataFrame
        updated client locations
    fac_df : geopandas.GeoDataFrame
        updated facility locations
    """
    
    fillers = [[cli_df, 'cli2fac'], [fac_df, 'fac2cli']]
    
    for df, attr in fillers:
        df[sol_col] = df['desc_var'].map(getattr(model, attr))
        df[sol_col].fillna('closed', inplace=True)
    
    if print_solution:
        selected = fac_df[fac_df[sol_col] != 'closed']
        for idx in selected.index:
            print('')
            print(selected.loc[idx, 'desc_var'], 'serving:',
                  selected.loc[idx, sol_col])
    
    return cli_df, fac_df


def service_area_object(fac_df, cli_df, f, style=None, smoother=10):
    """Create a geometric object for service area representation.
    either a `libpysal.cg.alpha_shape_auto()` object or a
    `shapely.geometry.MultiLineSting` object.
    
    Parameters
    ----------
    fac_df : geopandas.GeoDataFrame
        GeoDataFrame of facility locations.
    cli_df : geopandas.GeoDataFrame
        GeoDataFrame of client locations.
    f : str
        facility decision variable name.
    style : str or dict
        method for plotting service area. Default is None.
        If dict see `create_network_service`.
    smoother : float or int
        buffer (meters). Default is 10.
    
    Returns
    -------
    service_area : geopandas.GeoDataFrame
        contents are a geometric object representing
        facility service area
    """
    
    if style == 'spider_lines':
        web = []
        f_point = fac_df[fac_df.desc_var==f].geometry.squeeze()
        for c_point in cli_df.geometry:
            web.append(LineString((f_point, c_point)))
        
        service_area = [[MultiLineString(web)]]
    
    elif style == 'concave_hull':
        
        # client location coordinates
        c_array = np.array(cli_df.geometry.apply(lambda pt:
                         [pt.x, pt.y]).squeeze().tolist())
        
        # facility location coordinates
        f_array = np.array(fac_df[fac_df.desc_var==f].geometry.apply(\
                         lambda pt: [pt.x, pt.y]).squeeze())
        
        # coordinates of all location in the set
        pts_array = np.vstack((c_array, f_array))
        
        # create alpha shape (concave hull)
        ccv = cg.alpha_shape_auto(pts_array, step=4)
        service_area = [ccv.buffer(smoother)]
    
    elif type(style) == dict:
        model = style['model']
        net = style['net']
        service_area = create_network_service(model, net, f,
                                              cli_df, fac_df)
        service_area = [service_area]
    
    # instantiate geodataframe
    service_area = gpd.GeoDataFrame(service_area,
                                    columns=["geometry"])
    
    return service_area


def create_network_service(model, net, f, cgdf, fgdf):
    """network-based service area from service to clients
    
    Parameters
    ----------
    model : pulp.pulp.LpProblem
    net : spaghetti.network.Network
    f : str
        facility decision variable name.
    cgdf : geopandas.GeoDataFrame
        GeoDataFrame of client locations.
    fgdf : geopandas.GeoDataFrame
        GeoDataFrame of facility locations.
    
    Returns
    -------
    service_lines : list
        shapely.geometry.MultiLineString of
        network service area.
    """
    
    
    def create_service_lines(model, net):
        """
        Returns
        -------
        serv_lines : list
            shapely LineString objects representing network service area
        """
        
        serv_lines = []
        
        for (fac,cli), (cloc,floc) in model.cli2fac_sp.items():
            if fac != f:
                continue
            
            snap = 'snapped_coords'
            start = fgdf.loc[(fgdf.desc_var == fac), snap].squeeze()
            end = cgdf.loc[(cgdf.desc_var == cli), snap].squeeze()
            
            # if points snapped to different network segments
            if (cloc,floc) != ('snap', 'snap'):
                line = [Point(net.node_coords[floc])]
                for node_id in net.alldistances[cloc][1][floc]:
                    line.append(Point(net.node_coords[node_id]))
                line = line + [Point(net.node_coords[cloc])]
            
            # if points snapped to the same network segment
            else:
                line = []
            
            # add start and end points
            line.insert(0, start)
            line.append(end)
            serv_lines.append(LineString(line))
            
        # convert indivdual lines to one set of multilines
        serv_lines = [MultiLineString(serv_lines)]
        
        return serv_lines
    
    service_lines = create_service_lines(model, net)
    return service_lines


def distance_analytics(mdls, df_only=False, style_only=False,
                       save_as_latex=None, idx_name=None, color='green',
                       ltx_round=2):
    """create stylized dataframe visualization of distance analytics
    
    Parameters
    ----------
    mdls : list
        all modeling scenarios
    df_only : bool
        return only the analytics dataframe
    style_only : bool
        return only the stylized dataframe
    save_as_latex : str
        file path and name. Default is None.
    idx_name : str
        dataframe index name.
    color : str
        the background color gradient. Default is 'green'.
    ltx_round : int
        round to this many decimal places. Default is 3.
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        distance analytics matrix
    style : pandas.io.formats.style.Styler
        style dataframe view
    """
    
    model_names_base, model_name_map = _map_names(mdls)
    boiler =  ' to assigned facility'
    stats = {'abs_min': 'Absolute min dist'+boiler,
             'abs_max': 'Absolute max dist'+boiler,
             'mean_min': 'Mean min dist per client'+boiler,
             'mean_max': 'Mean max dist per client'+boiler,
             'mean_means': 'Mean of mean dists per client'+boiler,
             'mean_stds': 'Mean of StD dists per client'+boiler}
    
    # instantiate dataframe
    df = gpd.GeoDataFrame()
    df[idx_name] = model_names_base
    
    for n in list(stats.keys()):
        df[n] = np.nan
    
    # calculate stat for each model
    for m in mdls:
        ab_mins = []
        ab_maxs = []
        mean_min = []
        mean_max = []
        means = []
        stds = []
        
        for f, cs in m.fac2cli.items():
            row = np.array([m.fac2iloc[f]])
            cols = np.array([m.cli2iloc[c] for c in cs])
            
            try:
                dists = m.cij[row, cols[None, :]]
            except IndexError:
                continue
            
            # stats
            ab_mins.append(dists.min())
            ab_maxs.append(dists.max())
            mean_min.append(dists.mean())
            mean_max.append(dists.mean())
            means.append(dists.mean())
            stds.append(dists.std())
        
        # fill cells
        calcs = [np.array(ab_mins).min(),
                 np.array(ab_maxs).max(),
                 np.array(mean_min).min(),
                 np.array(mean_max).max(),
                 np.array(means).mean(),
                 np.array(stds).mean()]
        
        label_calc = {k: calcs[idx] for idx, k\
                                    in enumerate(list(stats.keys()))}
        
        for k, v in label_calc.items():
            df.loc[(df[idx_name] == model_name_map[m.name]), k] = v
    
    # set index as stat names
    df.set_index(idx_name, inplace=True)
    
    # stylize
    cm = sns.light_palette(color, as_cmap=True, reverse=True)
    style = df.style.set_caption(stats)\
                    .background_gradient(axis=0, cmap=cm)
    
    # convert df entry to latex
    if save_as_latex:
        ltx_df = df.copy()
        ltx_df = ltx_df.round(ltx_round)
        df_as_latex = _as_latex(ltx_df, color, idx_name=idx_name,
                                table_type='distance')
        
        # save out as printed version, not simple string
        with open('%s.tex' % save_as_latex, 'w') as tex:
            print(df_as_latex, file=tex)
    
    if df_only:
        return df
    
    if style_only:
        return style
    
    return df, style


def coverage_analytics(mdls, df_only=False, style_only=False,
                       save_as_latex=None, idx_name=None, color='green',
                       ltx_round=2):
    """create stylized dataframe visualization of coverage analytics
    
    Parameters
    ----------
    mdls : list
        all modeling scenarios
    df_only : bool
        return only the analytics dataframe
    style_only : bool
        return only the stylized dataframe
    save_as_latex : str
        file path and name. Default is None.
    idx_name : str
        dataframe index name.
    color : str
        the background color gradient. Default is 'green'.
    ltx_round : int
        round to this many decimal places. Default is 3.
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        distance analytics matrix
    style : pandas.io.formats.style.Styler
        style dataframe view
    """
    
    model_names_base, model_name_map = _map_names(mdls)
    stats = {'cli_not_cov': 'Client locations *not* covered'}
    
    most_coverage = 0
    for m in mdls:
        covered = max(m.ncov2ncli)
        if covered > most_coverage:
            most_coverage = covered
    coverage = list(range(1, most_coverage+1))
    
    for cv in coverage:
        stats['cov_by_%s' % cv] = 'Clients covered by %s facility' % cv
    
    # instantiate dataframe
    df = gpd.GeoDataFrame()
    df[idx_name] = model_names_base
    
    for n in list(stats.keys()):
        df[n] = np.nan
    
    # calculate stat for each model
    for m in mdls:
        calcs = [int(m.ncov2ncli[0])]
        
        for cv in coverage:
            try:
                c = m.ncov2ncli[cv]
            except KeyError:
                c = 0
            calcs.append(c)
        
        label_calc = {k: calcs[idx] for idx, k\
                                    in enumerate(list(stats.keys()))}
        
        for k, v in label_calc.items():
            df.loc[(df[idx_name] == model_name_map[m.name]), k] = v
    
    # set index as stat names
    df.set_index(idx_name, inplace=True)
    
    # stylize
    cm = sns.light_palette(color, as_cmap=True, reverse=True)
    style = df.style.set_caption(stats)\
                    .background_gradient(axis=0, cmap=cm)
    
    # convert df entry to latex
    if save_as_latex:
        ltx_df = df.copy()
        ltx_df = ltx_df.round(ltx_round)
        df_as_latex = _as_latex(ltx_df, color, idx_name=idx_name,
                                table_type='coverage')
        
        # save out as printed version, not simple string
        with open('%s.tex' % save_as_latex, 'w') as tex:
            print(df_as_latex, file=tex)
    
    if df_only:
        return df
    
    if style_only:
        return style
    
    return df, style


def solution_analytics(mdls, model_type=None, df_only=False, style_only=False,
                       save_as_latex=None, idx_name=None, good='green',
                       bad='red', ltx_round=2):
    """create stylized dataframe visualization of solution analytics
    
    Parameters
    ----------
    mdls : list
        all modeling scenarios
    model_type : str
        name of location model.
    df_only : bool
        return only the analytics dataframe
    style_only : bool
        return only the stylized dataframe
    save_as_latex : str
        file path and name. Default is None.
    idx_name : str
        dataframe index name.
    good : str
        the background color gradient and the Optimal solution
        label color. Default is 'green'.
    bad : str
        the Infeasible solution label color. Default is 'red'.
    ltx_round : int
        round to this many decimal places. Default is 3.
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        distance analytics matrix
    style : pandas.io.formats.style.Styler
        style dataframe view
    """
    
    
    def _highlight_sol(s, optimal, infeasible):
        """highlight solution optimality in pandas.DataFrame.
        """
        
        return ['background-color: %s' % optimal\
                if v == 'Optimal'\
                else 'background-color: %s' % infeasible\
                if v == 'Infeasible'\
                else '' for v in s]
    
    model_names_base, model_name_map = _map_names(mdls)
    stats = {'status': 'Solution status',
             'gap': 'Infeasibility gap',
             'sol_time': 'Solve time (minutes)'}
    
    if mdls[0].time_limit:
        stats['time_lim'] = 'Time limit for solution (minutes)'
    
    if model_type == 'pmp':
        stats['obj'] = 'Minimized population-weighted distance'
    
    if model_type == 'pcp':
        stats['obj'] = 'Minimized worst-case distance'
    
    if model_type == 'lscp':
        stats['obj'] = 'Facilities needed for objective (may be relaxed)'
    
    if model_type == 'mclp':
        stats['pop_not_cov'] = 'Total population *not* covered'
        stats['pct_not_cov'] = 'Percentage of population *not* covered'
    
    stats.update({
                  'constrs': 'Constraints count',
                  'vars': 'Variables count',
                  'dims': 'Cost matrix dimensions',})
    
    if model_type != 'lscp':
        stats['ndv=p'] = 'Selected facility variables equal p'
    
    # instantiate dataframe
    df = gpd.GeoDataFrame()
    df[idx_name] = model_names_base
    
    for n in list(stats.keys()):
        df[n] = np.nan
    
    # calculate stat for each model
    for m in mdls:
        
        # fill cells
        calcs = [m.status,
                 m.gap,
                 m.solve_time]
        
        if m.time_limit:
            calcs.append(m.time_limit)
        
        if model_type == 'pmp':
            calcs.append(m.obj_val)
        
        if model_type == 'pcp':
            calcs.append(m.obj_val)
        
        if model_type == 'lscp':
            calcs.append(int(m.obj_val))
        
        if model_type == 'mclp':
            calcs.append(m.ai_sum - m.obj_val)
            calcs.append(100. - m.perc_served)
        
        dims = (m.cij.shape[0], m.cij.shape[1])
        n_constrs, n_vars = element_calculator(model_type,
                                               dims[0],
                                               cols=dims[1])
        dims = str(dims)
        calcs.extend([n_constrs, n_vars, dims])
        
        if model_type != 'lscp':
            dv_match_p = True
            if len(m.fac2cli) != m.p:
                dv_match_p = False
            calcs.append(dv_match_p)
        
        label_calc = {k: calcs[idx] for idx, k\
                                    in enumerate(list(stats.keys()))}
        
        for k, v in label_calc.items():
            df.loc[(df[idx_name] == model_name_map[m.name]), k] = v
    
    # set index as stat names
    df.set_index(idx_name, inplace=True)
    
    # convert to integer
    int_stats = ['pop_not_cov', 'constrs', 'vars']
    for stat in int_stats:
        if stat in df.columns:
            df[stat] = df[stat].astype(int)
    if model_type == 'lscp':
         df['obj'] = df['obj'].astype(int)
    
    # stylize
    cm = sns.light_palette(good, as_cmap=True, reverse=True)
    style = df.style.set_caption(stats)\
                    .apply(_highlight_sol, args=(good, bad,))\
                    .background_gradient(axis=0, cmap=cm)
    
    # convert df entry to latex
    if save_as_latex:
        ltx_df = df.copy()
        ltx_df = ltx_df.round(ltx_round)
        df_as_latex = _as_latex(ltx_df, good, idx_name=idx_name,
                                bad=bad, table_type='solution')
        
        # save out as printed version, not simple string
        with open('%s.tex' % save_as_latex, 'w') as tex:
            print(df_as_latex, file=tex)
    
    if df_only:
        return df
    
    if style_only:
        return style
    
    return df, style


def selection_analytics(mdls, df_only=False, style_only=False, transpose=True,
                        save_as_latex=None, idx_name=None, full_color='green',
                        partial_color='orange', ltx_round=2):
    """create stylized dataframe visualization of
    selected decision variable analytics.
    
    Parameters
    ----------
    mdls : list
        all modeling scenarios
    df_only : bool
        return only the analytics dataframe. Default is False.
    style_only : bool
        return only the stylized dataframe. Default is False.
    transpose : bool
        transpose dataframe [True]. Default is True. When [False]
        a color gradient is added to the sum and percent columns.
    save_as_latex : str
        file path and name. Default is None.
    idx_name : str
        dataframe index name.
    full_color : str
        the full set color (selected facilities that serve clients).
        Default is 'green'.
    partial_color : str
        the partial set color (selected facilities that do not serve
        clients). Default is 'orange'.
    ltx_round : int
        round to this many decimal places. Default is 2.
    
    Returns
    -------
    df : geopandas.GeoDataFrame
        variable selection matrix
    style : pandas.io.formats.style.Styler
        style dataframe view
    """
    
    
    def _highlight_membership(s, full_color, partial_color):
        """highlight set membership in pandas.DataFrame.
        """
        
        if full_color == 'green':
            full_color_background = 'limegreen'
        if partial_color == 'orange':
            partial_color_background = 'orange'
        
        return ['background-color: %s' % full_color_background\
                if v == '$\in$'\
                else 'background-color: %s' % partial_color_background\
                if v == '$\\varnothing$'\
                else '' for v in s]
    
    
    model_names_base, model_name_map = _map_names(mdls)
    
    # set index and columns in empty dataframe
    var_index = [v.name for v in mdls[0].fac_vars]
    df = gpd.GeoDataFrame(index=var_index, columns=model_names_base)
    ncols = df.shape[1]
    df.index.name = idx_name
    
    # if site was selected in a model label with
    # latex symbol for 'element of a set' ($\in$)
    for m in mdls:
        for f in df.index:
            if f in list(m.fac2cli.keys()):
                if len(m.fac2cli[f]) == 0:
                    df.loc[f, model_name_map[m.name]] = '$\\varnothing$'
                else:
                    df.loc[f, model_name_map[m.name]] = '$\in$'
    
    # label all other cells with latex ($\\notin$)
    df.fillna('$\\notin$', inplace=True)
    in_syms = ['$\in$','$\\varnothing$']
    for idx in df.index:
        sel = df.loc[idx][df.loc[idx].isin(in_syms)].shape[0]
        df.loc[idx, '$\sum$'] = sel
        df.loc[idx, '$\%$'] = (float(sel) / float(ncols)) * 100.
    df['$\sum$'] = df['$\sum$'].astype(int)
    
    # transpose and set subset
    if transpose:
        df = df.T
        df.loc['$\%$', :] = df.loc['$\%$', :].astype(float).round(ltx_round)
        subset = None
    else:
        subset = ['$\sum$', '$\%$']
    
    # stylize
    cm = sns.light_palette(full_color, as_cmap=True)
    caption = 'Set membership by modeling scenario: '\
              + '[green]=selected facility serving clients; '\
              + '[orange]=selected facility serving no clients'
    style = df.style\
              .set_caption(caption)\
              .apply(_highlight_membership, args=(full_color, partial_color))\
              .background_gradient(cmap=cm, subset=subset)
    
    # convert df entry to latex
    if save_as_latex:
        ltx_df = df.copy()
        df_as_latex = _as_latex(ltx_df, full_color,
                                idx_name=idx_name, partial_color=partial_color,
                                table_type='selection')
        
        # save out as printed version, not simple string
        with open('%s.tex' % save_as_latex, 'w') as tex:
            print(df_as_latex, file=tex)
    
    if df_only:
        return df
    
    if style_only:
        return style
    
    return df, style


def _map_names(mdls):
    """internal function for mapping long to base model names
    
    Parameters
    ----------
    mdls : see `solution_analytics()`.
    
    Returns
    -------
    base : list
        base name for models to use as dataframe index
    name_map : dict
        {long: base} naming mapper.
    """
    
    long = [m.name for m in mdls]
    base = [m.name.split('_')[0] for m in mdls]
    name_map = dict(zip(long, base))
    
    return base, name_map


def _as_latex(df, full_color, idx_name, partial_color=None,
              bad=None, table_type=None):
    """convert dataframe to latex style table. This requires several
    steps and is a bit tedious.
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        dataframe of results
    full_color : see `selection_analytics()`.
    idx_name : see `selection_matrix()` or `analytics_matrix()`.
    partial_color : see `selection_analytics()`.
        Default is None.
    bad : see `solution_analytics()`.
    
    table_type : str
        {'distance', 'selection', 'coverage', 'solution'}.
        Default is None.
    
    Returns
    -------
    latex : str
        latex table style syntax
    """
    
    df.insert(0, idx_name, df.index)
    df.reset_index(drop=True, inplace=True)
    
    latex = df.to_latex(escape=False, index=False)
    latex = latex.replace('_', '\_')
    
    
    # set table full 'in' color
    if full_color == 'green':
        html_full_color = 'ABEBC6'
    
    if table_type == 'selection':
        latex = latex.replace('$\in$',
                              '\cellcolor[HTML]{%s}$\in$'\
                              % html_full_color)
    
    if table_type == 'solution':
        latex = latex.replace('Optimal',
                              '\cellcolor[HTML]{%s}Optimal'\
                              % html_full_color)
    
    # set table partial 'in' color
    if partial_color == 'orange':
        html_partial_color = 'FAD7A0'
    if table_type == 'selection' and partial_color:
        latex = latex.replace('$\\varnothing$',
                              '\cellcolor[HTML]{%s}$\\varnothing$'\
                              % html_partial_color)
    
    # set table bad solution color
    if bad == 'red':
        html_bad_color = 'F1948A'
    if table_type == 'solution' and bad:
        latex = latex.replace('Infeasible',
                              '\cellcolor[HTML]{%s}Infeasible'\
                              % html_bad_color)
    
    # set up table format
    if table_type == 'selection':
        headers = len(df.columns) - 2
        left, center = 'l' * headers, 'c' * (headers - 1)
        
        latex = latex.replace('\\begin{tabular}{%srr}' % left,
                              '\\begin{tabular}{l|%s|rr}' % center)
        
    if table_type in ['coverage', 'distance', 'solution']:
        headers = len(df.columns)
        right = 'r' * (headers - 1)
        
        latex = latex.replace('\\begin{tabular}{l%s}' % right,
                              '\\begin{tabular}{l|%s}' % right)
    
    return latex



###############################################################################
##########################    PLOTTING FUNCTIONS    ###########################
###############################################################################



def plotter(fig=None, base=None, plot_aux=None, buffered=None, model=None,
            pt1_size=None, pt2_size=None, plot_res=None, save_fig=None,
            title=None, area=None, census_geo=None, dv_col=None,
            w_serv=True, no_space=False, cls_plt=True, fdv='desc_var',
            sa_style='spider_lines', fs='x-large', figsize=(10,10)):
    """ Top-level scenario plotter for location analytics.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        complete figure to plot. Default is None.
    base : matplotlib.axes._subplots.AxesSubplot
        individual axis to plot. Default is None.
    plot_aux : dict
        model data parameters dataframes to plot keyed by
        descriptive names. Default is None.
    plot_res : dict
        model data results dataframes to plot keyed by
        descriptive names. Default is None.
    buffered : see
        buffer distance from roads segments in `plot_base`. Default is None.
    census_geo : str
        spatial unit of census geography.
    pt1_size : float or float
        size of points to plot. `pt1_size` should always be the
        larger between `pt2_size` and `pt1_size`. Default is None.
    pt2_size : float or float
        size of points to plot. Default is None.
    model : pulp.pulp.LpProblem
        location model
    area : str
        location of model
    title : str
        title. Default is None.
    sa_styleplot : {str, dict}
        method for plotting service area. Default is spider_lines.
        Option is concave_hull, or dict for network area.
    dv_col : str
        decision variable column name. Default is None.
    w_serv : bool
        add service area to legend [True]. Default is True.
    no_space : bool
        no space between legend entries [True]. Default is False.
    fdv : str
        facility label decision variable column. Default is 'desc_var'.
    cls_plt : bool
        close ax following plotcall [True]. Default is True.
    fs : str
        font size. Default is 'x-large'.
    figsize : tuple
        Figure size for plot. Default is (10,10).
    save_fig : str
        filename. Default is False.
    
    Returns
    -------
    add_to_legend : list
        items to add to legend
    """
    
    if area == 'Leon_FL':
        figsize = (25,20)
    
    for_multiplot = True
    if not fig and not base:
        for_multiplot = False
        fig, base = plt.subplots(1, 1, figsize=figsize)
        
        if area == 'Test_Tract_Leon_FL':
            base.set_xlim(620800, 625800)
    
    # add title
    if not for_multiplot:
        if title:
            title = ' '.join(title.split('_'))
        if model:
            if title:
                title += ' - ' + model.locmod
            else:
                title = model.name
        base.set_title(title, size=20)
    else:
        base.set_title(model.name, size=20)
    
    # plot non-results data
    if plot_aux:
        for k, df in plot_aux.items():
            if k == 'census_polys':
                df.plot(ax=base, color='k', alpha=.05, linestyle='--',
                        linewidth=2, zorder=1)
            if k == 'streets':
                if area == 'Test_Tract_Leon_FL':
                    street_width = 1
                else:
                    street_width = 2
                if type(sa_style) == dict:
                    street_width = .25
                df.plot(ax=base, lw=street_width, color='k', zorder=1)
            if k == 'buffer':
                df.plot(ax=base, color='y', lw=.25, alpha=.25, zorder=1)
            if k == 'cli_tru':
                if plot_res:
                    df = df[df[dv_col] == 'closed']
                    if df.empty:
                        continue
                    psize = pt2_size/6.
                    pcolor = 'k'
                else:
                    n_cli = df.shape[0]
                    psize = pt1_size
                    pcolor = 'b'
                df.plot(ax=base, markersize=psize,
                        edgecolor='k', color=pcolor)
            if k == 'fac_tru':
                if plot_res:
                    df = df[df[dv_col] == 'closed']
                    if df.empty:
                        continue
                    psize = pt2_size
                    pcolor = 'k'
                    pmarker = '*'
                else:
                    n_cli = df.shape[0]
                    psize = pt1_size
                    pcolor = 'r'
                    pmarker = 'o'
                df.plot(ax=base, markersize=psize*2,
                        edgecolor='k', color=pcolor,
                        marker=pmarker)
                n_fac = df.shape[0]
            if k == 'cli_snp':
                df.plot(ax=base, markersize=pt2_size,
                        edgecolor='k', color='b', alpha=.75)
            if k == 'fac_snp':
                df.plot(ax=base, markersize=pt2_size,
                        edgecolor='k', color='r', alpha=.75)
        add_to_legend = list(plot_aux.keys())
    else:
        add_to_legend = None
    
    # plot results data
    if plot_res:
        dv_colors = dv_colorset(plot_res['fac_var'][fdv])
        
        # facilities
        df = plot_res['fac_var'][plot_res['fac_var'][dv_col] != 'closed']
        alpha = .5
        
        # decision variable info for legend
        dvs_to_leg = {}
        
        # plot facilities
        for desc_var in df.desc_var:
            fac = df[df.desc_var == desc_var]
            fac.plot(ax=base, marker='*', markersize=pt1_size*3.,
                     alpha=.8, zorder=3, edgecolor='k',
                     color=dv_colors[desc_var])
            
            # update decision variable info with set color
            dvs_to_leg[desc_var] = {'color':dv_colors[desc_var]}
        
        # plot clients & service areas
        for f, c in model.fac2cli.items():
            fc = plot_res['cli_var'][plot_res['cli_var'].desc_var.isin(c)]
            fc.plot(ax=base, markersize=50, edgecolor='k',
                    color=dv_colors[f], alpha=alpha, zorder=2)
            
            # update decision variable info with set client counts
            dvs_to_leg[f].update({'clients': fc.shape[0]})
            
            # create service area object
            service_area = service_area_object(df, fc, f, style=sa_style)
            
            if type(sa_style) == str:
                # set appropriate alpha and line widths
                if sa_style == 'spider_lines':
                    records = len(service_area.geometry.squeeze())
                    if records < 10:
                        _alpha = .5
                    elif records >= 10 and records <= 20:
                        _alpha = .4
                    elif records >= 20 and records <= 30:
                        _alpha = .3
                    else:
                        _alpha = .15
                else:
                    _alpha = .2
                
                service_area.plot(ax=base, edgecolor='k', alpha=_alpha,
                                  lw=pt2_size/6., color=dv_colors[f], zorder=1)
            else:
                service_area.plot(ax=base, alpha=.2, zorder=1,
                                  color=dv_colors[f], linewidth=10)
    
    else:
        dvs_to_leg = None
    
    # create a shell class to represent
    # FacilityLocationModel if not present
    if not model:
        try:
            model = _ShellModel(plot_aux,
                                ['cli_tru', 'fac_tru', 'census_polys'])
        except (TypeError, KeyError):
            model = None
    
    # if FacilityLocationModel present add extra
    # attributes for patch creation
    else:
        try:
            key = 'census_polys'
            attr_name = 'n_' + key[:3]
            setattr(model, attr_name, plot_aux[key].shape[0])
        except KeyError:
            pass
    
    if not for_multiplot:
        
        # create legend patches
        patches = create_patches(model=model, for_multiplot=for_multiplot,
                                 pt1_size=pt1_size, pt2_size=pt2_size,
                                 buffered=buffered, legend_aux=add_to_legend,
                                 dvs_to_leg=dvs_to_leg, census_geo=census_geo,
                                 area=area, w_serv=w_serv, no_space=no_space)
        
        add_legend(patches, for_multiplot=for_multiplot, fs=fs)
    
    add_north_arrow(base, area=area)
    
    add_scale(base, area=area, for_multiplot=for_multiplot)
    
    if save_fig:
        plt.savefig(save_fig+'.png', bbox_inches='tight')
    
    # if for a multiplot explicityly return items to add to legend
    if for_multiplot:
        return add_to_legend
    
    if cls_plt:
        plt.close()


class _ShellModel:
    """object to mimic `model` when not present
    """
    
    def __init__(self, plot_aux, keys):
        
        for key in keys:
            attr_name = 'n_' + key[:3]
            try:
                setattr(self,
                        attr_name,
                        plot_aux[key].shape[0])
            except KeyError:
                pass


def multi_plotter(models, plot_aux=None, plot_res=None, select=None,
                  title=None, area=None, census_geo=None, net=None,
                  w_serv=True, no_space=False, sa_style='spider_lines',
                  fs='x-large', figsize=(14,14), shape=(2,2)):
    """plot multiple base axes as one figure
    
    Parameters
    ----------
    models : list
        solved model objects
    select : dict
        facility-to-selection count lookup.
    shape : tuple
        dimension for subplot array. Default is (2,2).s
    sa_style : str
        see plotter(). Default is 'spider_lines'.
    net : SpaghettiNetwork
        used for create of network service area
    census_geo : see plotter()
    plot_aux : see plotter()
    plot_res : see plotter()
    title : see plotter()
    area : see plotter()
    figsize : see plotter()
    w_serv : see plotter()
    fs : see plotter()
    no_space : see plotter()
    """
    
    pt1_size, pt2_size = 300, 60
    
    # convert list of models to array 
    mdls = np.array(models).reshape(shape)
    fig, axarr = plt.subplots(mdls.shape[0], mdls.shape[1],
                              figsize=figsize,
                                      sharex='col',
                                      sharey='row')
    
    # add super title to subplot array
    plt.suptitle(title, fontsize=30)
    fig.subplots_adjust(hspace=0.1, wspace=0.005, top=.925)
    
    # create each subplot
    for i in range(mdls.shape[0]):
        for j in range(mdls.shape[1]):
            if net:
                sa_style = {'model':mdls[i,j], 'net':net}
            add_to_legend = plotter(base=axarr[i,j],
                                    plot_aux=plot_aux,
                                    plot_res=plot_res,
                                    model=mdls[i,j],
                                    pt1_size=pt1_size,
                                    pt2_size=pt2_size,
                                    area=area, sa_style=sa_style)
            axarr[i,j].set_aspect('equal')
    
    add_to_legend = set(add_to_legend)
    
    # decision variable color set
    dv_colors = dv_colorset(plot_res['fac_var'].desc_var)
    dvs_to_leg = {f: dv_colors[f] for m in models\
                                  for f in m.fac2cli.keys()}
    
    # set ordered dict of {iloc:fac_var, color, x-selected}
    # *** models[0] can be any of the solved models
    dvs_to_leg = {models[0].fac2iloc[k]:(k,v, select[k])\
                  for k, v in dvs_to_leg.items()}
    dvs_to_leg = OrderedDict(sorted(dvs_to_leg.items()))
    
    try:
        _shell_model = _ShellModel(plot_aux, ['census_polys'])
    except (TypeError, KeyError):
        _shell_model = None
    
    # create legend patches
    patches = create_patches(model=_shell_model,
                             pt1_size=pt1_size, pt2_size=pt2_size,
                             legend_aux=add_to_legend, dvs_to_leg=dvs_to_leg,
                             census_geo=census_geo, area=area, w_serv=w_serv,
                             no_space=no_space, cls_plt=False,
                             for_multiplot=True)
    
    add_legend(patches, fs=fs, for_multiplot=True)


def add_north_arrow(base, area=None):
    """add a north arrow to an axes
    
    Parameters
    ----------
    base : see plotter()
    """
    
    if area == 'Phoenix_grid':
        x, y = 221200, 267200
    elif area == 'Leon_FL':
        x, y = 580000,180000
    elif area == 'Test_Grid_Leon_FL':
        x, y = 1, 1
    elif area == 'Test_Tract_Leon_FL':
        x, y = 621400, 162400
    
    try:
        x, y
    except UnboundLocalError:
        err_msg = 'Study area not defined internally.'
        raise UnboundLocalError(err_msg)
    
    arw = 'rarrow, pad=0.25'
    
    bbox_props = dict(boxstyle=arw, fc='w', ec='k', lw=2, alpha=.75)
    
    base.text(x, y, '      z    ', bbox=bbox_props,
             fontsize='large',fontweight='heavy',
             ha='center', va='center', rotation=90)


def add_scale(base, area=None, for_multiplot=False):
    """add a scale arrow to an axes
    
    Parameters
    ----------
    base : see plotter()
    area : see plotter()
    for_multiplot : see plotter()
    """
    
    if area == 'Test_Grid_Leon_FL':
        return
    
    # set scale anchor
    if area == 'Test_Tract_Leon_FL':
        offset = 250
    else:
        offset = 75
    
    x, y = base.get_xlim()[0]+offset, base.get_ylim()[0]+offset
    
    # set scale distance and units
    if area == 'Phoenix_grid':
        distance, units = .25, 'km'
    elif area == 'Leon_FL':
        distance, units = 10, 'km'
    elif area == 'Test_Tract_Leon_FL':
        distance, units = .75, 'km'
    try:
        distance, units
    except UnboundLocalError:
        err_msg = 'Study area not defined internally.'
        raise UnboundLocalError(err_msg)
    
    if for_multiplot:
        fontsize = 'small'
    else:
        fontsize = 'medium'
    
    scale_text = '|    ~%s%s~    |' % (distance, units)
    
    bbox_props = dict(boxstyle='round, pad=.5',
                      fc='w', ec='0.5', alpha=0.7)
    
    base.text(x, y, scale_text, fontsize=fontsize,
              fontstyle='italic', bbox=bbox_props)


def add_legend(patches, for_multiplot=False, fs=None):
    """Add a legend to a plot
    
    Parameters
    ----------
    patches : list
        legend handles matching plotted items
    for_multiplot : bool
        flag for multiplot. changes legend location
    fs : see plotter()
    """
    
    if for_multiplot:
        anchor = (1.1, 1.65)
    else:
        anchor = (1.005, 1.015)
    
    legend = plt.legend(handles=patches, loc='upper left',
                        fancybox=True, framealpha=.85,
                        bbox_to_anchor=anchor, fontsize=fs)
    legend.get_frame().set_facecolor('white')


def dv_colorset(dvs):
    """decision variables color set
    https://htmlcolorcodes.com/color-names/
    
    Parameters
    ---------
    dvs : geopandas.GeoSeries
        facility decision variables
    
    Returns
    -------
    dv_colors : dict
        decision variable to set color lookup
    """
    
    dv_colors = ['fuchsia', 'mediumseagreen', 'blueviolet',
                 'darkslategray', 'lightskyblue', 'saddlebrown',
                 'cyan', 'darkgoldenrod', 'limegreen', 'orange',
                 'coral', 'mediumvioletred', 'darkcyan',
                 'tomato', 'deeppink']
    
    dv_colors = {desc_var:dv_colors[idx] for idx, desc_var\
                 in enumerate(dvs)}
    
    return dv_colors


def create_patches(model=None, pt1_size=None, pt2_size=None, area=None,
                   buffered=None, legend_aux=None, dvs_to_leg=None,
                   for_multiplot=False, style=None, census_geo=None,
                   w_serv=True, no_space=False):
    """create all patches to add to the legend.
    
    Parameters
    ----------
    for_multiplot : bool
        for a single plot (True), or multiplot (False).
        Default is False.
    model : see plotter()
    pt1_size : see plotter()
    pt2_size : see plotter()
    buffered : see plotter()
    census_geo : see ploter()
    legend_aux : see plotter()
    dvs_to_leg : see plotter()
    area : see plotter()
    w_serv : see plotter()
    no_space : see plotter()
    
    Returns
    -------
    patches : list
        legend handles matching plotted items
    """
    
    if pt1_size:
        ms1 = float(pt1_size)/6.
    if pt2_size:
        ms2 = float(pt2_size)/8.
    
    # all patches to add to legend
    patches = gpd.GeoDataFrame()
    
    spacer = mpatches.Patch(fc='w', ec='w', linewidth=0,
                            alpha=.0, label='')
    
    # streets -- always plot
    if area == 'Test_Tract_Leon_FL':
        street_width = 1
    else:
        street_width = 2
    strs = mlines.Line2D([], [], color='k', label='Streets',
                         alpha=1, linewidth=street_width)
    patches = patches.append([['spacer', spacer],
                              ['strs', strs]])
    
    # non-results data
    if legend_aux:
        if 'buffer' in legend_aux:
            label = 'Street buffer (%sm)' % buffered
            strbuff = mpatches.Patch(fc='y', ec='k', linewidth=2,
                                     alpha=.5, label=label)
            patches = patches.append([['spacer', spacer],
                                      ['strbuff', strbuff]])
        
        if 'census_polys' in legend_aux:
            label = census_geo.title() + ' ($n$=%s)' % model.n_cen
            if census_geo == 'blocks':
                label = 'Pop. ' + label
            cenpoly = mpatches.Patch(fc='k', ec='k', alpha=.05,
                                     linestyle='--', linewidth=2,
                                     label=label)
            patches = patches.append([['spacer', spacer],
                                      ['cenpoly', cenpoly]])
        
        if 'census_cents' in legend_aux:######################################
            pass
        
        
        if 'cli_tru' in legend_aux:
            try:
                if dvs_to_leg:
                    pcolor = 'k'
                    msize = ms2/3.
                    plabel = 'Uncovered Clients '\
                             + '($n$=%s)' % model.n_cli_uncov
                else:
                    pcolor = 'b'
                    msize = ms1
                    plabel = 'Synth. Households ($n$=%s)' % model.n_cli
                cli_tru = mlines.Line2D([], [], color=pcolor,
                                        marker='o', ms=msize,
                                        linewidth=0, alpha=.8,
                                        markeredgecolor='k',
                                        label=plabel)
                patches = patches.append([['spacer', spacer],
                                          ['cli_tru', cli_tru]])
                
            except AttributeError:
                pass
        
        if 'fac_tru' in legend_aux:
            if dvs_to_leg:
                pcolor = 'k'
                msize = ms2
                pmarker = '*'
                no_fac = model.n_fac - len(list(model.fac2cli.keys()))
                plabel = 'Unselected Facilities ($n$=%s)' % no_fac
            
            else:
                pcolor = 'r'
                msize = ms1*2
                pmarker = 'o'
                plabel = 'Fire Stations '\
                         + '($n$=%s)' % model.n_fac
                if area == 'Test_Tract_Leon_FL':
                    plabel = 'Synth. ' + plabel
            
            fac_tru = mlines.Line2D([], [], color=pcolor,
                                    marker=pmarker, ms=msize,
                                    markeredgecolor='k',
                                    linewidth=0, alpha=1,
                                    label=plabel)
            patches = patches.append([['spacer', spacer],
                                      ['fac_tru', fac_tru]])
        
        
        if 'cli_snp' in legend_aux:
            label = 'Households snapped to network'
            cli_snp = mlines.Line2D([], [], color='b', marker='o',
                                    ms=ms2, linewidth=0, alpha=1,
                                    markeredgecolor='k', label=label)
            patches = patches.append([['spacer', spacer],
                                      ['cli_snp', cli_snp]])
            
            
        if 'fac_snp' in legend_aux:
            label = 'Fire Stations snapped to network'
            fac_snp = mlines.Line2D([], [], color='r', marker='o',
                                    ms=ms2, markeredgecolor='k',
                                    linewidth=0, alpha=1,
                                    label=label)
            patches = patches.append([['spacer', spacer],
                                      ['fac_snp', fac_snp]])
    
    patches = patches.append([['spacer', spacer]])
    
    # results data for single plot
    if dvs_to_leg and not for_multiplot:
        
        # add facility, client, and service area patches to legend
        for k, v in dvs_to_leg.items():
            fdv_label = 'Fire Station %s' % k
            if area == 'Test_Tract_Leon_FL':
                fdv_label = 'Synth. ' + fdv_label
            fdv = mlines.Line2D([], [], color=v['color'], marker='*',
                                ms=ms1/2., markeredgecolor='k',
                                linewidth=0, alpha=.8, label=fdv_label)
            cdv_label = 'Clients of %s ' % k \
                        + '($n$=%s)' % v['clients']
            cdv = mlines.Line2D([], [], color=v['color'], marker='o',
                                ms=ms1/6., markeredgecolor='k',
                                linewidth=0, alpha=.5, label=cdv_label)
            serv_label = '%s service area' % k
            serv = mpatches.Patch(fc=v['color'], ec=v['color'], linewidth=2,
                                  alpha=.25, label=serv_label)
            
            if w_serv:
                patches = patches.append([['spacer', spacer],
                                          ['fdv', fdv],
                                          ['cdv', cdv],
                                          ['serv', serv],
                                          ['spacer', spacer]])
            else:
                patches = patches.append([['spacer', spacer],
                                          ['fdv', fdv],
                                          ['cdv', cdv],
                                          ['spacer', spacer]])
    
    # results data for multiplot
    if dvs_to_leg and for_multiplot:
        for idx, (k, v, n) in dvs_to_leg.items():
            fdv = mlines.Line2D([], [], color=v, marker='*', ms=ms1/2,
                                markeredgecolor='k', linewidth=0,
                                alpha=.8, label='%s ($n$=%s)' % (k,n))
            patches = patches.append([['spacer', spacer],
                                      ['fdv', fdv],
                                      ['spacer', spacer]])
    
    patches = patches.values
    
    if no_space:
        patches = patches[patches[:,0] != 'spacer']
    
    return list(patches[:,1])


def image_stitcher(mdls, outfile, dims, plt_dir=''):
    """stitch together a mosaic of 4 .png files
    using PIL -- Pillow
    
    Parameters
    ----------
    mdls : list
        all modeling scenarios
    outfile : str
        file name to save
    dims : tuple
        (rows, columns) mosaic dimensions
    plt_dir : str
        plots directory. Default is ''.
    
    Returns
    -------
    mosaic : PIL.Image.Image
        mosaiced image
    """
    
    rows, cols = dims[0], dims[1]
    
    # read in all images
    imgs = [Image.open('%s%s.png' % (plt_dir, m.name)) for m in mdls]
    
    # verify image size equality
    img_sizes = [img.size for img in imgs]
    if not all(img_sizes):
        raise RunTimeError('Images do not have equal dimensions. '\
                           + 'Stitch manually.')
    
    # extract width and height information
    width, height = img_sizes[0][0], img_sizes[0][1]
    
    # create a new image instance
    mosaic = Image.new('RGB', (width*rows, height*cols))
    
    # set location anchors
    locs = []
    for r in range(rows):
        for c in range(cols):
            locs.append([r*width, c*height])
    
    # iterate over images and set each image to its anchor
    for idx, img in enumerate(imgs):
        mosaic.paste(img, locs[idx])
    
    # save out
    mosaic.save('%s%s.png' % (plt_dir, outfile))
    
    return mosaic
