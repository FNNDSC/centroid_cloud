#!/usr/bin/env python2.7

import  random
import  numpy   as np
import  getopt
import  sys
import  math

_numPoints          = 100
_f_xoffset          = 0.0
_f_yoffset          = 0.0
_f_xscale           = 1.0
_f_yscale           = 1.0
_str_distribFunc    = 'None'
_str_outFile        = 'cloud.txt'

Gstr_synopsis= """

    NAME
    
        cloud_create.py
        
    SYNOPSIS
    
        cloud_create.py    [-s <outFile>]            \\
                           [-n <points>]             \\
                           [-O <xoffset,yoffset>]    \\
                           [-S <xscale,yscale>]      \\
                           [-d <distribFunc>]        
                           
    DESCRIPTION
    
        'cloud_create.py' creates a "cloud" of random points. By default,
        values in X and Y are spread between 0.0 and 1.0.
        
        Several additional parameters control the position, extend, and
        distribution of the cloud. The '-O' offset controls where the
        cloud's lower left boundary is, and the '-S' scale controls the
        scaling on the random range. Finally, the '-d' controls the 
        shape of the cloud.
        
    ARGUMENTS
    
        -s <outFile>
        The output file to contain the cloud.
    
        -n <points>
        The number of points in the cloud. Default is %d.
    
        -O <xoffset,yoffset>
        The lower left corner of the cloud. By default this is the
        coordinate space origin, (0; 0). If set, must specify both
        the xoffset and the yoffset.
        
        -S <xscale,yscale>
        The extent of the cloud. Defaults to (1, 1);
        
        -d <distribFunc>
        The function to control the shape of the cloud. Default is 'None', 
        i.e. no function.
        
        Valid functions include:
        
        o 'scaleX':
            Each y-value of the cloud is scaled by the corresponding
            x-value. This creates a right-angled triangle.
        o 'linearX': 
            Each y-value of the cloud is increased by the current   
            x-value. This creates a cloud that is skewed "upwards".
        o 'linearXscaled':
            Similar to 'linearX', but each value is also scaled by the
            current x-value. this creates a linear funnel.
        o 'quadX':   
            Each y-value of the cloud is increased by the square
            of the current x-value. This creates a quadratically
            skewed cloud.
        o 'quadXscaled':   
            Similar to 'quadX', but each value is also scaled by the
            current x-value. This creates a quadratic funnel.
        o 'circle': 
            If specified, the 'x' dimension is taken to be angle, and scaled
            by 2*pi; and the 'xscale' is not applicable. The 'y' dimension
            is taken to be radius, and scaled by 'yscale'.
        o 'shape':
            An internal shape (e.g. a cross) that can be used to test
            symmetry and asymmetry.

""" % (_numPoints)

_dict_distrib = {
    'None':             lambda p: np.array(p),
    'shape':            lambda p: np.array(p),
    'scaleX':           lambda p: np.column_stack((p[:,0],p[:,1] * p[:,0])),
    'linearX':          lambda p: np.column_stack((p[:,0],p[:,1] + p[:,0])),
    'linearXscaled':    lambda p: np.column_stack((p[:,0],(p[:,1] + p[:,0]) * p[:,0])),
    'quadX':            lambda p: np.column_stack((p[:,0],p[:,1] + p[:,0]**2)), 
    'quadXscaled':      lambda p: np.column_stack((p[:,0],(p[:,1] + p[:,0]**2) * p[:,0])),
    'circle':           lambda p: np.column_stack((p[:,1]*np.cos(p[:,0]),p[:,1]*np.sin(p[:,0])))
                 }


def synopsis_show():
    print "USAGE: %s" % Gstr_synopsis

try:
    opts, remargs   = getopt.getopt(sys.argv[1:], 'n:s:O:S:d:xh')
except getopt.GetoptError:
    sys.exit(1)

for o, a in opts:
    if (o == '-x' or o == '-h'):
        synopsis_show()
        sys.exit(1)
    if (o == '-s'):
        _str_outFile = a
    if (o == '-n'):
        _numPoints = int(a)
    if (o == '-O'):
        _str_XYoffset       = a
        _l_XYoffset         = _str_XYoffset.split(',')
        _f_xoffset          = float(_l_XYoffset[0])
        _f_yoffset          = float(_l_XYoffset[1])
    if (o == '-S'):
        _str_XYscale        = a
        _l_XYscale          = _str_XYscale.split(',')
        _f_xscale           = float(_l_XYscale[0])
        _f_yscale           = float(_l_XYscale[1])
    if (o == '-d'):
        _str_distribFunc    = a

def internals_print():
    print   "-n numPoints   = %d"  % _numPoints
    print   "-O offset      = %f, %f"  % (_f_xoffset, _f_yoffset)
    print   "-S scale       = %f, %f"  % (_f_xscale, _f_yscale)
    print   "-d distribFunc = %s"  % _str_distribFunc
    

_M_cloud    = np.random.random( (_numPoints, 2) )
_M_cloud   *= np.array( (_f_xscale, _f_yscale) )
if _str_distribFunc != 'circle':
    _M_cloud   += np.array( (_f_xoffset, _f_yoffset) )
    if _str_distribFunc == 'shape':
        _M_left     = _M_cloud.copy() + np.array( (0+_f_xoffset, -1+_f_yoffset) )
        _M_center   = _M_cloud.copy() + np.array( (1+_f_xoffset,  0+_f_yoffset) )
        _M_right    = _M_cloud.copy() + np.array( (2+_f_xoffset,  0+_f_yoffset) )
        _M_top      = _M_cloud.copy() + np.array( (1+_f_xoffset,  1+_f_yoffset) )
        _M_bottom   = _M_cloud.copy() + np.array( (-1+_f_xoffset, -2+_f_yoffset) )
        _M_cloud    = np.vstack((_M_left, _M_center, _M_right, _M_top, _M_bottom))
#        _M_cloud    = np.vstack((_M_left, _M_center, _M_right, _M_top))
else:
    _M_cloud[:,0] *= 2*math.pi
   
_M_cloud    = _dict_distrib[_str_distribFunc](_M_cloud)   
if _str_distribFunc == 'circle':
    _M_cloud   += np.array( (_f_xoffset, _f_yoffset) )
    
np.savetxt(_str_outFile, _M_cloud, fmt='%10.5f', delimiter='\t', newline='\n')

        
