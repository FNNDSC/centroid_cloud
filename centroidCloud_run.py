#!/usr/bin/env python2.7

from    pylab   import *
import  numpy   as np
import  getopt
import  sys

from    C_centroidCloud import *

Gstr_synopsis = """
    
    NAME
 
        centroidCloud_run.py
        
    SYNOPSIS
    
        centroidCloud_run.py    -m <cloudFile>       \\
                                [-r <rotations>      \\
                                -s <stdWidth>        \\
                                -z <zorder>          \\
                                -n                   \\
                                -a                   \\
                                -A                   \\
                                -e -d                \\
                                -x -h]            
                                
    DESCRIPTION
    
        centroidCloud_run.py is a simple "driver" for the C_centroidCloud
        class. This class accepts a file defining a cloud of points in a 
        2D space, and then returns a polygon of points defining the projection
        of the cloud deviation along a set of rotated axes.
        
        In the default case, there are 90 rotations corresponding to a 1 degree
        sweep through the 1st Cartesian quadrant.
        
    ARGUMENTS
    
        -m <cloudFile>
        The cloud file to read.

        -r <rotations>
        The number of rotations, evenly spread out between 0...90 degrees.

        -s <stdWidth>
        The width of the standard deviation cloud

        -z <zorder>
        Change the 'z' ordering. By default, the confidence boundary is drawn
        "on top" of the point plots with a default zorder=3. By specifying
        a different zorder (such as '-z 1'), the boundary can be drawn below
        the point plots.
        
        -a
        If specified, will turn ON axis equal flag.

        -A
        If specified, will turn ON asymmetrical deviations flag.

        -n 
        If specified, will turn OFF clould normalization.

        Default is to always normalize. Without normalizing, 
        rotational skew can occur. Normalization adds extra operations 
        (and hence time) to the projection calculations. This extra time
        is minimal, and normalization should probably be used in all cases.

        -a
        If specified, will set axis ranges equal. In cases where deviation is
        much more along one axis than another, this will result in a very
        "thin" boundary polygon.
        
        -e
        If specified, print an extent report that gives for each rotation
        the X, Y, XY, and X+Y projection extent. This is an approximation
        for the max/min extent angle for the given cloud.
        
        -d
        If specified, pass a debug flag to the C_centroidCloud class that
        triggers additional reporting.
                        
        -x or -h
        Print this help page.
        
"""

Gstr_cloudMatrix            = 'cloud.txt'
Gb_saveFig                  = False
Gb_debugCloud               = False
Gb_extentReport             = False
Gb_normalize                = True
Gb_asymmetricalDeviations   = False
Gb_axisEqual                = False
Gi_zorder                   = 3
G_numRotations              = 90
G_f_stdWidth                = 0.5

def synopsis_show():
    print "USAGE: %s" % Gstr_synopsis

def deviation_plot(al_points, str_fillColor = 'red', str_edgeColor = 'black'):
    poly    = Polygon(al_points,
                        facecolor = str_fillColor,
                        edgecolor = str_edgeColor, zorder=Gi_zorder)
    gca().add_patch(poly)
    return poly

try:
    opts, remargs   = getopt.getopt(sys.argv[1:], 'hxm:s:r:dez:naA')
except getopt.GetoptError:
    sys.exit(1)

for o, a in opts:
    if (o == '-x' or o == '-h'):
        synopsis_show()
        sys.exit(1)
    if (o == '-m'):
        Gstr_cloudMatrix            = a
    if (o == '-s'):
        G_f_stdWidth                = float(a)
    if (o == '-r'):
        G_numRotations              = int(a)
    if (o == '-z'):
        Gi_zorder                   = int(a)
    if (o == '-n'):
        Gb_normalize                = False
    if (o == '-a'):
        Gb_axisEqual                = True
    if (o == '-A'):
        Gb_asymmetricalDeviations   = True
    if (o == '-d'):
        Gb_debugCloud               = True
    if (o == '-e'):
        Gb_extentReport             = True

C_cloud         = C_centroidCloud(file='%s' % Gstr_cloudMatrix, 
                                  stdWidth  = G_f_stdWidth,
                                  rotations = G_numRotations)
C_cloud.normalize(Gb_normalize)
C_cloud.asymmetricalDeviations(Gb_asymmetricalDeviations)
C_cloud.debug(Gb_debugCloud)
C_cloud.confidenceBoundary_find()
M_cloud         = C_cloud.cloud()
l_polygonPoints = C_cloud.boundary()

figure()
if Gb_axisEqual:
    axis('equal')
grid() 

p1,     = plot(M_cloud[:,0], M_cloud[:,1], color='#FF0000', marker='*', ls='None',
               zorder = 1)
poly    = deviation_plot(l_polygonPoints, 'green') 
 
if Gb_extentReport: print C_cloud.projectionExtent_report()
show()

