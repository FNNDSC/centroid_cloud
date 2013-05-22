#!/usr/bin/env python

from    pylab       import  *
import  numpy       as      np
import  sys
from   _common      import systemMisc as misc
import  math

class C_centroidCloud:
    """
    Determines the confidence polygon boundary of a cloud of points
    by projecting extent along a rotating set of axes.
    
    Assumes 2D clouds, although should be expandable to N-D without
    fundamental changes. Eg, 3D is a series of 2D polygons calculated
    on a set rotating planes.
    """
    
    
    # 
    # Class member variables -- if declared here are shared
    # across all instances of this class
    #
    mdictErr = {
        'Keys'          : {
            'action'        : 'initializing base class, ',
            'error'         : 'it seems that no member keys are defined.',
            'exitCode'      : 10},
        'noCloud'          : {
            'action'        : 'initializing base class, ',
            'error'         : 'it seems that no cloud data was provided.',
            'exitCode'      : 11},
        'Save'              : {
            'action'        : 'attempting to pickle save self, ',
            'error'         : 'a PickleError occured',
            'exitCode'      : 12},
        'SaveMat'           : {
            'action'        : 'attempting to save MatLAB friendly spectrum, ',
            'error'         : 'an IOerror occured',
            'exitCode'      : 13},
        'Load'              : {
            'action'        : 'attempting to pickle load object, ',
            'error'         : 'a PickleError occured',
            'exitCode'      : 14},
        'Rotations'         : {
            'action'        : 'intializing base class, ',
            'error'         : 'numRotations must be at least 1.',
            'exitCode'      : 20},
    }

    def dprint( self, level, str_txt ):
        """
            Simple "debug" print... based on verbosity level.
        """
        if level <= self.m_verbosity: print(str_txt)

    def verbosity_set( self, level ):
        self.m_verbosity = level

    def error_exit( self,
                            astr_key,
                            ab_exitToOs=1
                            ):
        print("FATAL ERROR")
        print("\tSorry, some error seems to have occurred in <%s>::%s" \
                        % ( self.__name__, self.__proc__ ))
        print("\tWhile %s" % C_centroidCloud.mdictErr[astr_key]['action'])
        print("\t%s" % C_centroidCloud.mdictErr[astr_key]['error'])
        print("")
        if ab_exitToOs:
            print("Returning to system with error code %d" % \
                            C_centroidCloud.mdictErr[astr_key]['exitCode'])
            sys.exit( C_centroidCloud.mdictErr[astr_key]['exitCode'] )
        return C_centroidCloud.mdictErr[astr_key]['exitCode']

    def fatal( self, astr_key, astr_extraMsg="" ):
        if len( astr_extraMsg ): print(astr_extraMsg)
        self.error_exit( astr_key )

    def warn( self, astr_key, astr_extraMsg="" ):
        b_exitToOS = 0
        if len( astr_extraMsg ): print(astr_extraMsg)
        self.error_exit( astr_key, b_exitToOS )

    def std2mean(self, aM, av_mean):
        """
        Calculates deviation of <aM> about to <av_mean>, which might
        be different than the actual mean of <aM>.
        """
        dims    = len(av_mean)
        v_dev   = np.copy(av_mean)
        N       = aM.shape[0]
        for dim in np.arange(0, dims):
            f_dev       = sqrt(np.sum(1.0/float(N) * (aM[:,dim] - av_mean[dim])**2))
            v_dev[dim]  = f_dev
        return v_dev

    def stats_calc(self, aM, adict_stats):
        """
        Assumes that the cloud <aM> is in row order, i.e.
             x1    y1
             x2    y2
               ...
             xn     yn
        
        Dimensionality of space is equal to number of cols

        Places stats (per column) in adict_stats.
        
        Returns the stats dictionary.
        """
        adict_stats['mean']     = np.mean(aM, 0)
        adict_stats['std']      = np.std(aM, 0)
        adict_stats['stdpos']   = np.std(aM, 0)
        adict_stats['stdneg']   = np.std(aM, 0)
        adict_stats['min']      = np.min(aM, 0)
        adict_stats['max']      = np.max(aM, 0)
        adict_stats['range']    = np.max(aM, 0) - np.min(aM,0)

        # Calculate the stdpos and stdneg extent
        dims = len(adict_stats['mean'])
        for dim in np.arange(0, dims):
            # pos
            if(self._str_centerMean == 'subset'):
                v_p = np.std(aM[aM[:,dim]>=adict_stats['mean'][dim]],0)
            if(self._str_centerMean == 'original'):
                v_p = self.std2mean(aM[aM[:,dim]>=adict_stats['mean'][dim]], 
                          adict_stats['mean'])
            adict_stats['stdpos'][dim] = v_p[dim]
            # neg
            if(self._str_centerMean == 'subset'):
                v_n = np.std(aM[aM[:,dim]<adict_stats['mean'][dim]],0)
            if(self._str_centerMean == 'original'):
                v_n = self.std2mean(aM[aM[:,dim]<adict_stats['mean'][dim]], 
                          adict_stats['mean'])
            adict_stats['stdneg'][dim] = v_n[dim]

        return adict_stats
                
    def rot_2D(self, aM, af_theta):
        """
        Return the 2D rotation of a cloud of 2D points, aM, 
        about angle af_theta (rad):
        
        aM is assumed to be of form:

           +-        -+
           | x1    y1 |
           | x2    y2 |
           | x3    y3 |
           |    ...   | 
           | xn    yn |
           +-        -+
           
        and the rotation is calculated according to:
           
        +-                               -+ +-              -+
        | cos(af_theta)    -sin(af_theta) | | x1  x2 ... xn  |
        | sin(af_theta)     cos(af_theta) | | y1  y2 ... yn  |
        +-                               -+ +-              -+
       
        i.e. 2x2 x 2xN with return 2xN (transposed):

           +-          -+
           | x1r    y1r |
           | x2r    y2r |
           | x3r    y3r |
           |    ...     | 
           | xnr    ynr |
           +-          -+
           
        
        """
        M_rot       = np.zeros( (2,2) )
        M_rot[0,0]  = math.cos(af_theta);   M_rot[0,1] = -math.sin(af_theta)
        M_rot[1,0]  = math.sin(af_theta);   M_rot[1,1] =  math.cos(af_theta)
        M_Ctr       = aM.transpose()
        M_Crot      = np.dot(M_rot, M_Ctr)
        return M_Crot.transpose()

    def cloudSpace_normalize(self, aM_cloud, o_adict_stats):
        '''
        This method "normalizes" the input cloud to a unit space along its
        basis dimensions. This reduces rotational skew in the original space
        which is particularly apparent when the original space deviations
        beteen different dimensions are very large.

        The normalized cloud is returned.
        '''
        self.stats_calc(aM_cloud, o_adict_stats)
        _M              = (aM_cloud - o_adict_stats['min']) / \
                            (o_adict_stats['max'] - o_adict_stats['min'])
        return _M


    def cloudSpace_denormalize(self, aM_cloud, adict_stats):
        '''
        This method "de-normalizes" the input cloud, based on
        the values passed in adict_stats which holds the stats
        of the original space.

        Note, this method overwrites the passed arguments!
        '''
        _M = aM_cloud * \
                    (adict_stats['max'] - adict_stats['min']) + \
                    adict_stats['min']
        return _M
    
    def projectionsOnAxes_find(self, adict_stats):
        """
        Determines the 2*nD points of the projections along the coordinate axes,
        in the frame of coordinate axes.

        The asymmetricalDeviations flag controls per-dimensional pos/neg 
        deviation calculations.
        
        Returns: ((x_min, x_max), (y_min, y_max), ... )
        """
        dims = len(adict_stats['mean'])
        v_ret = np.zeros( (dims), dtype='object')
        for dim in np.arange(0, dims):
            if self._b_asymmetricalDeviations:
                v_projection = np.array([adict_stats['mean'][dim] - 
                                         adict_stats['stdneg'][dim]*self._f_dev,
                                         adict_stats['mean'][dim] + 
                                         adict_stats['stdpos'][dim]*self._f_dev])
            else:
                v_projection = np.array([adict_stats['mean'][dim] - 
                                         adict_stats['std'][dim]*self._f_dev,
                                         adict_stats['mean'][dim] + 
                                         adict_stats['std'][dim]*self._f_dev])
            v_ret[dim] = v_projection
        return v_ret
    
    def projectionsOnAxes_rotate(self, av_projection, af_theta, a_d_stats):
        """
        Rotate projections in a reference axis frame by <af_theta>.
        
        Returns a rotated matrix:
            +-                           -+            \
            | dim0_min_projection_rotated |    <-- QIII \  
            | dim0_max_projection_rotated |    <-- Q1     > 2D
            | dim1_min_projection_rotated |    <-- QIV  /
            | dim1_max_projection_rotated |    <-- QII /
            |            ...              |           / 
            +-                           -+
        """
        
        # First, pack the line endpoints into a coordinate-pair matrix
        # structure
        M_p             = np.zeros( (2*self._M_Cdimensionality, 
                                     self._M_Cdimensionality) )
        index           = 0
        for dim in np.arange(0, self._M_Cdimensionality):
            for endpoint in [0, 1]:
                M_p[index, dim] = av_projection[dim][endpoint]
                index += 1
        M_center        = np.tile(a_d_stats['mean'],
                                  (2* self._M_Cdimensionality, 1))
        M_mask          = np.logical_not(M_p).astype(float)
        M_centerMask    = M_center * M_mask
        M_p             = M_p + M_centerMask
        # Now rotate the end points by <af_theta>
        M_p_rot = self.rot_2D(M_p, af_theta)
        return M_p_rot

    def projectionExtent_find(self, av_projections):
        """
        Simply calculate "distance" between the min/max points on 
        a projection for all the dimensions.
        
        Return the extent in an array, indexed by dimension
        """
        v_extent    = np.zeros(self._M_Cdimensionality)
        for dim in np.arange(0, self._M_Cdimensionality):
            v_extent[dim] = av_projections[dim][1] - av_projections[dim][0]
        return np.absolute(v_extent)
    
    def extent_init(self, v_extent):
        self._dict_extent2D['X']['min']     = v_extent[0] 
        self._dict_extent2D['X']['max']     = v_extent[0] 
        self._dict_extent2D['Y']['min']     = v_extent[1] 
        self._dict_extent2D['Y']['max']     = v_extent[1] 
        self._dict_extent2D['XY']['min']    = v_extent[0] * v_extent[1]
        self._dict_extent2D['XY']['max']    = v_extent[0] * v_extent[1]
        self._dict_extent2D['X+Y']['min']   = v_extent[0] + v_extent[1]
        self._dict_extent2D['X+Y']['max']   = v_extent[0] + v_extent[1]

    def extent_process(self, key, f_extent, rotation):
        if self._dict_extent2D[key]['min'] >= f_extent:
            self._dict_extent2D[key]['min'] = f_extent
            self._dict_extent2D[key]['minAngle'] = rotation
        if self._dict_extent2D[key]['max'] <= f_extent:
            self._dict_extent2D[key]['max'] = f_extent
            self._dict_extent2D[key]['maxAngle'] = rotation
              
    def projectionExtent_report(self):
        """
        Return a string that reports on the min/max extents on the X, Y axes
        """
        str_ret = ""
        for key in self._dict_extent2D.keys():
            for str_val in ['min', 'max']: 
                str_ret += "%s %3s projection @ angle = %10.5f @ %5.2f\n" % \
                            (str_val, key, 
                             self._dict_extent2D[key][str_val],
                             self._dict_extent2D[key]['%sAngle' % str_val])
        return str_ret

    def confidenceBoundary_find(self):
        """
        The main controlling function for finding the confidence boundary
        of a cloud.
        
        For each rotation angle in the internal dictionary, rotate the cloud,
        find the projections on the xy axes, then rotate the projections back.
        
        Store a counterclock-wise set of boundary polygon points:
        Q1->Q2->Q3->Q4 in the self._l_boundary list
        
        NOTE:
        * Only properly debugged for planar (i.e. 2D) boundaries.
        """
        self.stats_calc(self._M_C, self._d_origCloudStats)
        if self._b_normalizeCloudSpace:
            self._M_C = self.cloudSpace_normalize(self._M_C, self._d_origCloudStats)
        for rotation in self._rotationKeys:
            f_angle = self._dict_rotationVal[rotation]
            f_rad   = np.deg2rad(f_angle)
            # First, rotation the cloud by -r_rad:
            neg_C   = self.rot_2D(self._M_C, -f_rad)
            # Now find the projections on the standard x/y axis:
            self.stats_calc(neg_C, self._d_rotatedCloudStats[rotation])
            self._d_rotatedCloudStats[rotation]['rotation'] = rotation
            v_projections = self.projectionsOnAxes_find(self._d_rotatedCloudStats[rotation])
            v_extent = self.projectionExtent_find(v_projections)
            if self._b_normalizeCloudSpace:
                v_extent = self.cloudSpace_denormalize(v_extent, self._d_origCloudStats)
                v_extent = np.absolute(v_extent)
            if self._b_debug:
                print("%f %f %f %f %f" % (self._dict_rotationVal[rotation],
                                        v_extent[0], v_extent[1], 
                                        v_extent[0] + v_extent[1],
                                        v_extent[0] * v_extent[1]))
            if rotation == 0: self.extent_init(v_extent)
            for key in self._dict_extent2D.keys():
                if key == 'X':      self.extent_process('X', v_extent[0], 
                                        self._dict_rotationVal[rotation])
                if key == 'Y':      self.extent_process('Y', v_extent[1], 
                                        self._dict_rotationVal[rotation])
                if key == 'X+Y':    self.extent_process('X+Y', 
                                        v_extent[0] + v_extent[1], 
                                        self._dict_rotationVal[rotation])
                if key == 'XY':     self.extent_process('XY', 
                                        v_extent[0] * v_extent[1], 
                                        self._dict_rotationVal[rotation])
            M_p = self.projectionsOnAxes_rotate(v_projections, f_rad, \
                                                self._d_rotatedCloudStats[rotation])
            if self._b_normalizeCloudSpace:
                M_p = self.cloudSpace_denormalize(M_p, self._d_origCloudStats)
            self._dict_Q[1][rotation] = M_p[1,:]
            self._dict_Q[2][rotation] = M_p[3,:]
            self._dict_Q[3][rotation] = M_p[0,:]
            self._dict_Q[4][rotation] = M_p[2,:]
        for quadrant in range(1, 5):
            for rotation in self._rotationKeys:
                self._l_boundary.append(self._dict_Q[quadrant][rotation])
        if self._b_normalizeCloudSpace:
            self._M_C   = self.cloudSpace_denormalize(self._M_C, self._d_origCloudStats)
        #self.minmaxComponent_analyze()

        
    def initialize(self):            
        """
        (Re)builds internal dictionary structures. Typically called
        once the number of rotations has been set and/or a new cloud
        has been read.
        """
        if self._numRotations < 1:
            self.fatal('Rotations')
        
        self._M_Crows, self._M_Ccols = self._M_C.shape
        self._M_Cdimensionality = self._M_Ccols
        # The numRotations defines the number of rotations between 
        # 0 and 90 degrees (i.e. the first quadrant). This is used
        # to create a dictionary of rotationKeys and rotationVals
        self._rotationKeys      = range(0, self._numRotations)
        v_keys                  = np.array(self._rotationKeys)
        v_vals                  = v_keys * 90.0/self._numRotations

        for rotation in self._rotationKeys:
            self._d_rotatedCloudStats[rotation] = self._dict_stats.copy()
        self._dict_rotationVal  = misc.dict_init(self._rotationKeys, 
                                                 list(v_vals))
        self._dict_projection   = misc.dict_init(self._rotationKeys, 
                                        np.zeros( (self._M_Ccols, self._M_Ccols),
                                                  dtype = 'object'))
        # 2D Planar quadrant projections
        for quadrant in range(1,5):
            self._dict_Q[quadrant] = misc.dict_init(self._rotationKeys,
                                        np.zeros( (1,2) ))

    def dimensionality(self):
        '''
        Return the dimensionality of the cloud
        '''
        return self._M_Cdimensionality

    def boundary(self, *args):
        """
        Get/set the boundary points.
        """
        if len(args):
            self._l_boundary = args[0]
        else:
            return self._l_boundary

    def cloud(self, *args):
        """
        Get/set the cloud points.
        
        If set, trigger a re-initialization of system.
        """
        if len(args):
            self._M_C = args[0]
            self.initialize()
        else:
            return self._M_C
        
    def deviationWidth(self, *args):    
        """
        Get/set the deviation width.
        
        If set, trigger a re-initialization of system.
        """
        if len(args):
            self._f_dev = args[0]
            self.initialize()
        else:
            return self._f_dev

    def normalize(self, *args):
        """
        Get/set the normalization flag.
        """
        if len(args):
            self._b_normalizeCloudSpace = args[0]
        else:
            return self._b_normalizeCloudSpace

    def asymmetricalDeviations(self, *args):
        """
        Get/set the asymmetricalDeviations flag.
        """
        if len(args):
            self._b_asymmetricalDeviations = args[0]
        else:
            return self._b_asymmetricalDeviations
        
    def centerMean(self, *args):
        """
        Get/set the centerMean value.
        """
        if len(args):
            self._str_centerMean = args[0]
        else:
            return self._str_centerMean

    def debug(self, *args):    
        """
        Get/set the debugging flag.
        """
        if len(args):
            self._b_debug = args[0]
        else:
            return self._b_debug
    
    def rotations(self, *args):    
        """
        Get/set the number of rotations.
        
        If set, trigger a re-initialization of system.
        """
        if len(args):
            self._numRotations = args[0]
            self.initialize()
        else:
            return self._numRotations
       
    def __init__( self, *args, **kwargs ):
        self.__name__                   = 'C_centroidCloud'
        self.__proc__                   = "__init__"
        self._b_debug                   = False
        
        self._v_centroid                = np.zeros( (1, 2) )

        self._M_C                       = None
        self._M_Crows                   = -1
        self._M_Ccols                   = -1
        self._M_Cdimensionality         = 0
        
        self._dict_stats                = {
                'mean':         [],
                'std':          [],
                'stdpos':       [],
                'stdneg':       [],
                'min':          [],
                'max':          [],
                'range':        [],
                'rotation':     0.0
        }
        self._l_min                     = []
        self._l_max                     = []
        self._d_origCloudStats          = self._dict_stats.copy()
        self._d_rotatedCloudStats       = {}
        self._f_dev                     = 1.0
        self._b_asymmetricalDeviations  = False
        self._str_centerMean            = 'original'
        self._b_normalizeCloudSpace     = True

        self._str_file                  = ''
        self._numRotations              = 90
        self._b_deg                     = True
        for key, value in kwargs.iteritems():
            if key == 'file':
                self._str_file      = value
                self._M_C           = np.genfromtxt(self._str_file)
            if key == 'cloud':      self._M_C           = value
            if key == 'rotations':  self._numRotations  = value
            if key == 'stdWidth':   self._f_dev         = value
            if key == 'normalize':  self._b_normalizeCloudSpace = value

        if not self._M_C.any():
            self.fatal('noCloud')

        # The rotations are set/controlled by _numRotations which defines
        # the number of rotations between 0 and 90 degrees. The rotation values
        # are stored in _dict_rotationVal, indexed by _rotationKeys. In the
        # trivial (default) case, the rotation vals and keys are identical.
        self._rotationKeys          = []
        self._dict_rotationVal      = {}
        self._dict_projection       = {}
        
        # Planar polygon
        # The boundary points of the projections along the rotated axes and
        # within a given plane are stored in four dictionaries, one for each 
        # cartesian quadrant:
        #
        #              |
        #          Q2  |  Q1
        #     ---------+---------
        #          Q3  |  Q4
        #              |
        #
        self._dict_Q                = {1: {}, 2: {}, 3:{}, 4:{}}
        self._l_boundary            = []
        
        #
        #
        # Properties of the distribution on a single planar surface
        #
        self._dict_minMax2D         = {'min': 0.0, 'minAngle': 0.0, 
                                       'max': 0.0, 'maxAngle': 0.0}
        self._dict_extent2D         = {'X':     self._dict_minMax2D.copy(),
                                       'Y':     self._dict_minMax2D.copy(),
                                       'XY':    self._dict_minMax2D.copy(),
                                       'X+Y':   self._dict_minMax2D.copy()}
        
        
        self.initialize()
        
        
                                    
