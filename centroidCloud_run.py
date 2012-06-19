#!/usr/bin/env python

from    pylab import *
import  getopt

from    C_centroidCloud import *

Gstr_matrixType = 'rand.txt'
Gb_saveFig      = False

def synopsis_show():
    print "USAGE:"

def deviation_plot(M_P, str_fillColor = 'red', str_edgeColor = 'black'):
    f_meanX = np.mean(M_P[:,0])
    f_stdX  = np.std(M_P[:,0])
    f_meanY = np.mean(M_P[:,1])
    f_stdY  = np.std(M_P[:,1])
    rect    = Rectangle([f_meanX - f_stdX/2, f_meanY - f_stdY/2], 
                        f_stdX, f_stdY,
                        facecolor = str_fillColor,
                        edgecolor = str_edgeColor)
    gca().add_patch(rect)
    return rect

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

C_cloud         = C_centroidCloud(file='%s' % Gstr_matrixType)
C_cloud.confidenceBoundary_find()


