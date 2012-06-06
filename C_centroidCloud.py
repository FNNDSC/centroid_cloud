#!/usr/bin/env python

Gstr_matrixType = 'lower'
Gb_saveFig      = False

from    pylab import *
import  numpy as np
import  getopt
import  sys


class C_centroidCloud:
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

    def stats_determine(self):
        self._dict_stats['meanX']   = np.mean(self._M_C[:,0])
        self._dict_stats['meanY']   = np.mean(self._M_C[:,1])
        self._dict_stats['stdX']    = np.std(self._M_C[:,0])
        self._dict_stats['stdY']    = np.std(self._M_C[:,1])

    def stats_calc(self):
        # Assumes that the cloud is in row order, i.e.
        #    x1    y1
        #    x2    y2
        #      ...
        #    xn     yn
        #
        rows, cols = self._M_C.shape
        for dim in np.arange(0, cols):
            self._dict_stats['mean'] = np.mean(self._M_C, 0)
            self._dict_stats['std']  = np.std(self._M_C, 0)
            

    def __init__( self, *args, **kwargs ):
        self.__name__       = 'C_centroidCloud'
        self._v_centroid    = np.zeros( (1, 2) )
        self._angleSteps    = 90
        
        self._str_file      = ''
        for key, value in kwargs.iteritems():
            if key == 'file':   self._str_file = value

        self._M_C   = np.genfromtxt(self._str_file)
        
        self._dict_stats            = {}
        self._dict_stats['mean']    = []
        self._dict_stats['std']     = []
        self._f_dev                 = 1.0
        self._numRotations          = 90                            
