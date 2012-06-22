#!/usr/bin/env python2.7

from    pylab   import *
import  numpy   as np
import  getopt
import  sys

from    C_centroidCloud import *

Gstr_matrixType = 'rand.txt'
Gb_saveFig      = False

def synopsis_show():
    print "USAGE:"

def deviation_plot(al_points, str_fillColor = 'red', str_edgeColor = 'black'):
    poly    = Polygon(al_points,
                        facecolor = str_fillColor,
                        edgecolor = str_edgeColor)
    gca().add_patch(poly)
    return poly

try:
    opts, remargs   = getopt.getopt(sys.argv[1:], 'hxm:s')
except getopt.GetoptError:
    sys.exit(1)

for o, a in opts:
    if (o == '-x' or o == '-h'):
        synopsis_show()
        sys.exit(1)
    if (o == '-m'):
        Gstr_matrixType = a
    if (o == '-s'):
        Gb_saveFig = True

C_cloud         = C_centroidCloud(file='%s' % Gstr_matrixType, 
                                  stdWidth  = 0.5,
                                  rotations = 90)
C_cloud.confidenceBoundary_find()

M_cloud         = C_cloud.cloud()
l_polygonPoints = C_cloud.boundary()

figure()
grid() 

p1, = plot(M_cloud[:,0], M_cloud[:,1], color='#FF0000', marker='*', ls='None')
poly    = deviation_plot(l_polygonPoints) 
 
show()
