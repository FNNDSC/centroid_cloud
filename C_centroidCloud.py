#!/usr/bin/env python

from    pylab       import  *
import  numpy       as      np
import  sys
import  systemMisc  as      misc
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
            'exitCode'      : 14}
    }

    def dprint( self, level, str_txt ):
        """
            Simple "debug" print... based on verbosity level.
        """
        if level <= self.m_verbosity: print str_txt

    def verbosity_set( self, level ):
        self.m_verbosity = level

    def error_exit( self,
                            astr_key,
                            ab_exitToOs=1
                            ):
        print "%s:: FATAL ERROR" % self.mstr_obj
        print "\tSorry, some error seems to have occurred in <%s::%s>" \
                        % ( self.__name__, self.mstr_def )
        print "\tWhile %s" % C_centroidCloud.mdictErr[astr_key]['action']
        print "\t%s" % C_centroidCloud.mdictErr[astr_key]['error']
        print ""
        if ab_exitToOs:
            print "Returning to system with error code %d" % \
                            C_centroidCloud.mdictErr[astr_key]['exitCode']
            sys.exit( C_centroidCloud.mdictErr[astr_key]['exitCode'] )
        return C_centroidCloud.mdictErr[astr_key]['exitCode']

    def fatal( self, astr_key, astr_extraMsg="" ):
        if len( astr_extraMsg ): print astr_extraMsg
        self.error_exit( astr_key )

    def warn( self, astr_key, astr_extraMsg="" ):
        b_exitToOS = 0
        if len( astr_extraMsg ): print astr_extraMsg
        self.error_exit( astr_key, b_exitToOS )

    def stats_calc(self, aM):
        """
        Assumes that the cloud <aM> is in row order, i.e.
             x1    y1
             x2    y2
               ...
             xn     yn
        
        Dimensionality of space is equal to number of cols
        
        Returns
        """
        _dict_stats         = {}
        _dict_stats['mean'] = np.mean(aM, 0)
        _dict_stats['std']  = np.std(aM, 0)
        return _dict_stats
                
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
        M_Crot      = M_rot * M_Ctr
        return M_Crot.transpose()
    
    def projectionsOnAxes_find(self, adict_stats):
        """
        Determines the 2*nD points of the projections along the coordinate axes,
        in the frame of coordinate axes.
        
        Returns: ((x_min, x_max), (y_min, y_max), ... )
        """
        dims = len(adict_stats['mean'])
        v_ret = np.zeros(1, dims, dtype='object')
        for dim in np.arange(0, dims):
            v_projection = np.array([adict_stats['mean'][dim] - 
                                     adict_stats['std'][dim] * self._f_dev,
                                     adict_stats['mean'][dim] + 
                                     adict_stats['std'][dim] * self._f_dev])
            v_ret[dim] = v_projection # copy?
        return v_ret
    
    def projectionsOnAxes_rotate(self, av_projection, af_theta):
        """
        Rotate projections in a reference axis frame by <af_theta>.
        
        Returns a matrix:
            +-                           -+            \
            | dim0_min_projection_rotated |    <-- QIII \  
            | dim0_max_projection_rotated |    <-- Q1     > 2D
            | dim1_min_projection_rotated |    <-- QIV  /
            | dim1_max_projection_rotated |    <-- QII /
            |            ...              |           
            +-                           -+
        """
        
        # First, pack the line endpoints into a coordinate-pair matrix
        # structure
        dimensionality = len(av_projection)
        M_p             = np.zeros( (2*dimensionality, dimensionality) )
        index           = 0
        for dim in np.arange(0, dimensionality):
            for endpoint in [0, 1]:
                M_p[index, dim] = av_projection[dim][0, endpoint]
                index += 1
        # Now rotate the end points by <af_theta>
        M_p_rot = self.rot_2D(M_p, af_theta)
        return M_p_rot

    def confidenceBoundary_find(self):
        """
        The main controlling function for finding the confidence boundary
        of a cloud.
        
        For each rotation angle in the internal dictionary, rotate the cloud,
        find the projections on the xy axes, then rotate the projections back.
        
        Return a clock-wise set of boundary polygon points.
        """
        for rotation in self._rotationKeys:
            f_angle = self._dict_rotationVal[rotation]
            f_rad   = f_angle * np.pi / 180
            # First, rotation the cloud by -r_rad:
            neg_C   = self.rot_2D(self._M_C, -f_rad)
            # Now find the projections on the standard x/y axis:
            _dict_stats = self.stats_calc(neg_C)
            v_projections = self.projectionsOnAxes_find(_dict_stats)
            M_p = self.projectionsOnAxes_rotate(v_projections, f_rad)
            self._dict_Q1[rotation] = M_p[1,:]
            self._dict_Q2[rotation] = M_p[3,:]
            self._dict_Q3[rotation] = M_p[0,:]
            self._dict_Q4[rotation] = M_p[2,:]
    
    def initialize(self):            
        """
        (Re)builds internal dictionary structures. Typically called
        once the number of rotations has been set and/or a new cloud
        has been read.
        """
        self._M_Crows, self._M_Ccols = self._M_C.shape
        self._rotationKeys      = range(0, self._numRotations)
        self._dict_rotationVal  = misc.dict_init(self._rotationKeys, 
                                                 self._rotationKeys)
        self._dict_projection   = misc.dict_init(self._rotationKeys, 
                                        np.zeros( (self._M_Ccols, self._M_Ccols),
                                                  dtype = 'object'))
        # 2D Planar quadrant projections
        self._dict_Q1           = misc.dict_init(self._rotationKeys,
                                        np.zeros( (1, 2) ))
        self._dict_Q2           = misc.dict_init(self._rotationKeys,
                                        np.zeros( (1, 2) ))
        self._dict_Q3           = misc.dict_init(self._rotationKeys,
                                        np.zeros( (1, 2) ))
        self._dict_Q4           = misc.dict_init(self._rotationKeys,
                                        np.zeros( (1, 2) ))
        
    def __init__( self, *args, **kwargs ):
        self.__name__       = 'C_centroidCloud'
        self._v_centroid    = np.zeros( (1, 2) )
        
        self._str_file      = ''
        for key, value in kwargs.iteritems():
            if key == 'file':   self._str_file = value

        self._M_C                   = np.genfromtxt(self._str_file)
        self._M_Crows               = -1
        self._M_Ccols               = -1
        
        self._dict_stats            = {}
        self._dict_stats['mean']    = []
        self._dict_stats['std']     = []
        self._f_dev                 = 1.0

        # The rotations are set/controlled by _numRotations which defines
        # the number of rotations between 0 and 90 degrees. The rotation values
        # are stored in _dict_rotationVals, indexed by _rotationKeys. In the
        # trivial (default) case, the rotation vals and keys are identical.
        self._b_deg                 = True
        self._numRotations          = 90
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
        self._dict_Q1               = {}
        self._dict_Q2               = {}
        self._dict_Q3               = {}
        self._dict_Q4               = {}
        
        self.initialize()
        
        
                                    
