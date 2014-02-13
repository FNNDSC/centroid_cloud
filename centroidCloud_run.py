#!/usr/bin/env python2.7
from    __future__              import print_function

from    pylab                   import  *
import  itertools
import  scipy.stats             as      stats
import  numpy                   as      np
import  getopt
import  sys
import	os

from    C_centroidCloud         import  *
from    _common                 import systemMisc as misc

Gstr_synopsis = """
    
    NAME
 
        centroidCloud_run.py
        
    SYNOPSIS
    
        centroidCloud_run.py    -m <cloudFileLst>                       \\
                                [-r <rotations>                         \\
                                -s <stdWidth>                           \\
                                -z <zorder>                             \\
                                -n                                      \\
                                -a                                      \\
                                [ -A <centerMean> | -p <percentile>]    \\
                                -C <cloudColorLst> -K <kernelColorLst>  \\
                                -N <depth> -S <scale>                   \\
                                -e -d                                   \\
                                -x -h]            
                                
    DESCRIPTION
    
        centroidCloud_run.py is a simple "driver" for the C_centroidCloud
        class. This class accepts a file defining a cloud of points in a 
        2D space, and then returns a polygon of points defining the projection
        of the cloud deviation along a set of rotated axes.
        
        In the default case, there are 90 rotations corresponding to a 1 degree
        sweep through the 1st Cartesian quadrant.
        
    ARGUMENTS
    
        -m <cloudFileLst>
        A comma separated list of cloud files to read. Each cloud file will 
        be plotted and its centroid cloud determined.

        -r <rotations>
        The number of rotations, evenly spread out between 0...90 degrees.

        -s <stdWidth>
        The width of the standard deviation cloud

        -z <zorder>
        Change the 'z' ordering. By default, the confidence boundary is drawn
        "on top" of the point plots with a default zorder=3. By specifying
        a different zorder (such as '-z 1'), the boundary can be drawn below
        the point plots.
        
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

        -A <centerMean>
        If specified, will turn ON asymmetrical deviations flag. The
        <centerMean> controls the concept of "center" for the asymmetrical
        deviations. The following options are understood:
        
            'original' (default): std is calculated on a subset of original
                                  observations relative to the original mean
                                  of the main cloud.
                                  
            'subset':             std is calculated relative to the mean
                                  of the subset.
                                  
        -p <percentile>
        If specified, trigger a descriptive statistical analysis, using the
        passed <percentile> as upper and lower deviation from the center 
        mean.    

        -C <cloudColorLst> -K <kernelColorLst>
        A set of comma separated color specifiers, defining the colors of the
        point cloud and the kernel region respectively.
        
        -N <depth> -S <scale>  
        If specified, will calculate and plot a range of shifted clouds. This
        is only applicable for cases where more than one cloud has been
        specified as input. The -N <depth> denotes the number of square
        rotations to perform, and the -S <scale> the scale for each 
        radial rotation.

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

Gstr_cloudMatrixLst         = 'cloud.txt'
Gb_saveFig                  = False
Gb_debugCloud               = False
Gb_extentReport             = False
Gb_normalize                = True
Gb_asymmetricalDeviations   = False
Gb_usePercentiles           = False
G_f_percentile              = 25
Gstr_centerMean             = 'original'
Gb_axisEqual                = False
Gi_zorder                   = 3
G_numRotations              = 90
G_f_stdWidth                = 0.5
Gstr_cloudColorLst          = 'red,green,yellow'
Gstr_kernelColorLst         = 'blue,cyan,magenta'
Gb_rotateClouds             = False
G_depth                     = 3
G_f_depthScaleX             = 0.1
G_f_depthScaleY             = 0.1


def synopsis_show():
    print("USAGE: %s" % Gstr_synopsis)

def deviation_plot(al_points, str_fillColor = 'red', str_edgeColor = 'black'):
    poly    = Polygon(al_points,
                        facecolor = str_fillColor,
                        edgecolor = str_edgeColor, zorder=Gi_zorder)
    gca().add_patch(poly)
    return poly

try:
    opts, remargs   = getopt.getopt(sys.argv[1:], 'hxm:s:r:dez:naA:p:C:K:N:S:')
except getopt.GetoptError:
    sys.exit(1)

for o, a in opts:
    if (o == '-x' or o == '-h'):
        synopsis_show()
        sys.exit(1)
    if (o == '-m'):
        Gstr_cloudMatrixLst         = a
    if (o == '-C'):
        Gstr_cloudColorLst          = a
    if (o == '-K'):
        Gstr_kernelColorLst         = a
    if (o == '-N'):
        Gb_rotateClouds             = True
        G_depth                     = int(a)
    if (o == '-S'):
        G_f_depthScaleX             = a.split(',')[0]
        G_f_depthScaleY             = a.split(',')[1]
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
    if (o == '-p'):
        Gb_usePercentiles           = True
        G_f_percentile              = float(a)
    if (o == '-A'):
        Gb_asymmetricalDeviations   = True
        Gstr_centerMean             = a
        if(Gstr_centerMean != 'original' and Gstr_centerMean != 'subset'):
            Gstr_centerMean = 'original'
    if (o == '-d'):
        Gb_debugCloud               = True
    if (o == '-e'):
        Gb_extentReport             = True

v_relrotScale       = np.array((G_f_depthScaleX, G_f_depthScaleY))
l_cloudFile         = Gstr_cloudMatrixLst.split(',')
l_cloudColor        = Gstr_cloudColorLst.split(',')
l_kernelColor       = Gstr_kernelColorLst.split(',')
lC_cloud            = []
lM_cloud            = []
ll_polygonPoints    = []
p                   = []
poly                = []


# First read in clouds:
for cloud in range(0, len(l_cloudFile)):
    print("Reading %s" % l_cloudFile[cloud])
    lC_cloud.append(
                  C_centroidCloud(file      ='%s' % l_cloudFile[cloud], 
                                  stdWidth  = G_f_stdWidth,
                                  rotations = G_numRotations)        
                  )
    print("Appending cloud matrix...")
    lM_cloud.append(lC_cloud[cloud].cloud())
    print("Creating boundary list holder...")
    ll_polygonPoints.append([])
print("Preprocessing complete\n")    

# Now process the clouds, with possible iterative rotations
if Gb_rotateClouds:
    l_rotaryPoints  = misc.neighbours_findFast(2, G_depth)
else:   
    l_rotaryPoints  = [ np.array((0,0)) ]

l_rotaryPoints = l_rotaryPoints * v_relrotScale


for reltran in l_rotaryPoints:
    print("Relative shape translation: ", end="")
    print(reltran)
    figure()
    if Gb_axisEqual:
        axis('equal')
    grid() 
    for cloud in range(0, len(l_cloudFile)):
        if cloud:
            # all clouds except the base are displaced by the 
            # current translation
            print("translating cloud %s" % l_cloudFile[cloud])
            lC_cloud[cloud].cloud(lM_cloud[cloud] + reltran)
        print("Processing %s" % l_cloudFile[cloud])
        print("Normalizing...")
        lC_cloud[cloud].normalize(Gb_normalize)
        lC_cloud[cloud].usePercentiles(Gb_usePercentiles)
        lC_cloud[cloud].percentile(G_f_percentile)
        lC_cloud[cloud].asymmetricalDeviations(Gb_asymmetricalDeviations)
        lC_cloud[cloud].centerMean(Gstr_centerMean)
        print("Finding confidence boundary...")
        lC_cloud[cloud].confidenceBoundary_find()
        print("Extracting polygonPoints...")
        # ll_polygonPoints.append(lC_cloud[cloud].boundary())
        ll_polygonPoints[cloud] = lC_cloud[cloud].boundary()
        print("Plotting...")
        p.append(plot(lM_cloud[cloud][:,0], lM_cloud[cloud][:,1], 
                        color=l_cloudColor[cloud], marker='*', ls='None',
                        zorder = 1))
        poly.append(deviation_plot(ll_polygonPoints[cloud], l_kernelColor[cloud])) 
        if Gb_extentReport: print(C_cloud[cloud].projectionExtent_report())

    if len(l_cloudFile) > 1:
        l_combinations = list(itertools.combinations(range(len(l_cloudFile)), 2))
        print(l_combinations)
        for combination in l_combinations:
            g1  = combination[0]
            g2  = combination[1]
            _str_title = '%s-%s-x%f-y%f' % (
            		os.path.splitext(os.path.basename(l_cloudFile[g1])),
            		os.path.splitext(os.path.basename(l_cloudFile[g2])),
            		reltran[0], reltran[1]
            		)
            print("group 1 = %d" % g1)
            print("group 2 = %d" % g2)
            v_tstat, v_pval = stats.ttest_ind(lM_cloud[g1], lM_cloud[g2], equal_var = True)
            print(v_pval)
            f_pval  = np.amin(v_pval)
            print("p-val for comparison between group %d and %d is %f" % (g1, g2, f_pval))
    show()

