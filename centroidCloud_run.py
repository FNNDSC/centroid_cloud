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

from    shapely.geometry        import  Point           as sgPoint
from    shapely.geometry        import  Polygon         as sgPolygon
from    shapely.geometry        import  MultiPolygon    as sgMultiPolygon

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
                                -O <offsetX,offsetY>                    \\
                                -e -d -X                                \\
                                -x -h]            
                                
    DESCRIPTION
    
        centroidCloud_run.py is a simple "driver" for the C_centroidCloud
        class. This class accepts a file defining a cloud of points in a 
        2D space, and then returns a polygon of points defining the projection
        of the cloud deviation along a set of rotated axes. In the default 
        case, there are 90 rotations corresponding to a 1 degree sweep 
        through the 1st Cartesian quadrant.

        Although this driver was initially conceived as a single cloud
        processor, it has been extensively expanded to function also as
        multiple cloud analysis tool. If multiple clouds are passed as
        input, this program will plot all clouds and their confidence
        boundaries, as well as analyze all the p-val separations between
        all combinations of inputs.

        The script will also accept offsets to non-base clouds and will
        shift non-base clouds accordingly. Finally, it can also loop non-base
        clouds over square spirals about the base-cloud. The base-cloud is
        simply the first cloud input into the system.
        
    ARGUMENTS
    
        -m <cloudFileLst>
        A comma separated list of cloud files to read. Each cloud file will 
        be plotted and its centroid cloud determined.

        -r <rotations>
        The number of rotations for the confidence boundary for each cloud, 
        evenly spread out between 0...90 degrees.

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

        -X 
        If specified, plot convex hull boundaries.

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

        -O <offsetX,offsetY>                   
        If specified, will shift non-base clouds by offset. This is useful for
        once-off tests.

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
Gb_convexHullUse            = False
G_f_percentile              = 25
Gstr_centerMean             = 'original'
Gb_axisEqual                = False
Gi_zorder                   = 3
G_numRotations              = 90
G_f_stdWidth                = 0.5
Gstr_cloudColorLst          = 'red,green,yellow'
Gstr_kernelColorLst         = 'blue,cyan,magenta'
Gb_rotateClouds             = False
G_depth                     = 0
G_f_depthScaleX             = 0.1
G_f_depthScaleY             = 0.1
Gb_nonBaseOffset            = False
G_f_nonBaseOffsetX          = 0.0
G_f_nonBaseOffsetY          = 0.0

def synopsis_show():
    print("USAGE: %s" % Gstr_synopsis)

def deviation_plot(al_points, str_fillColor = 'red', str_edgeColor = 'black'):
    poly    = Polygon(al_points,
                        facecolor = str_fillColor,
                        edgecolor = str_edgeColor, zorder=Gi_zorder,
                        alpha     = 0.5)
    gca().add_patch(poly)
    return poly

def aspectRatio_square(aspect = 1):
    extent                    = np.array( (0, 6.20, 4.80, 0))
    print(extent)
    axes().set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def cloud_normalize(aM):
    '''
    For a cloud of N-dimensional points in <aM>, return a vector of normals
    to each point. This has the effect of reducing a high dimensional space
    to a single dimensional vector.
    '''
    _v_norm = []
    for p in aM:
        _v_norm.append(np.linalg.norm(p))
    return _v_norm

def convexHull_boundaryFind(ar_boundary):
    '''
    For a given np array <ar_boundary>, deterime the convex hull 
    (implicitly assuming 2D spaces).
    
    Basically, this builds a polygon, finds the convex hull, and
    translates back to an np.array.
    '''
    return np.asarray( sgPolygon(ar_boundary).convex_hull.exterior )

try:
    opts, remargs   = getopt.getopt(sys.argv[1:], 'hxm:s:r:dez:naA:p:C:K:N:S:O:X')
except getopt.GetoptError:
    sys.exit(1)

for o, a in opts:
    if (o == '-x' or o == '-h'):
        synopsis_show()
        sys.exit(1)
    if (o == '-X'):
        Gb_convexHullUse            = True
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
        G_f_depthScaleX             = float(a.split(',')[0])
        G_f_depthScaleY             = float(a.split(',')[1])
    if (o == '-O'):
        G_f_nonBaseOffsetX          = float(a.split(',')[0])
        G_f_nonBaseOffsetY          = float(a.split(',')[1])
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

v_reltranScale      = np.array((G_f_depthScaleX, G_f_depthScaleY))
l_cloudFile         = Gstr_cloudMatrixLst.split(',')
l_cloudColor        = Gstr_cloudColorLst.split(',')
l_kernelColor       = Gstr_kernelColorLst.split(',')
lC_cloud            = []
lM_cloudOrig        = []
lv_normOrig         = []
lM_cloud            = []
lv_norm             = []
ll_polygonPoints    = []
p                   = []
poly                = []


# First read in clouds:
for cloud in range(0, len(l_cloudFile)):
    print("Reading '%s'..." % l_cloudFile[cloud])
    lC_cloud.append(
                  C_centroidCloud(file      ='%s' % l_cloudFile[cloud], 
                                  stdWidth  = G_f_stdWidth,
                                  rotations = G_numRotations)        
                  )
    print("\tappending cloud matrix and norm vectors...")
    lM_cloudOrig.append(lC_cloud[cloud].cloud())
    lM_cloud.append(lC_cloud[cloud].cloud())
    lv_normOrig.append(cloud_normalize(lM_cloud[cloud]))
    lv_norm.append(lv_normOrig[cloud])
    print("\tcreating boundary list holder...")
    ll_polygonPoints.append([])
print("\nPreprocessing complete.\n")    

# Now process the clouds, with possible iterative rotations
if Gb_rotateClouds:
    l_rotaryPoints  = misc.neighbours_findFast(2, G_depth, includeOrigin=True)
else:   
    l_rotaryPoints  = [ np.array((0,0)) ]

v_nonBaseOffset = np.array( (G_f_nonBaseOffsetX, G_f_nonBaseOffsetY))
l_rotaryPoints  = l_rotaryPoints * v_reltranScale + v_nonBaseOffset

v_indoffset     = np.array( (G_depth, G_depth ) ) * 0.0;
M_pvalInd       = l_rotaryPoints.copy() + v_indoffset
v_pval          = np.zeros((len(M_pvalInd),1))
M_pvalflt       = np.column_stack( (M_pvalInd, v_pval) )
M_pval          = np.zeros((2*G_depth+1, 2*G_depth+1))
M_pvalmNorm     = np.zeros((2*G_depth+1, 2*G_depth+1))
M_pvalvNorm     = np.zeros((2*G_depth+1, 2*G_depth+1))

indexTotal      = 0
b_baseProcessed = False
for reltran in l_rotaryPoints:
    print(M_pval)
    figure(indexTotal)
    axis('auto')
    if Gb_axisEqual:
        axis('equal')
    grid() 
    for cloud in range(0, len(l_cloudFile)):
        print("\nProcessing '%s' containing %d points..." % \
                (l_cloudFile[cloud], len(lM_cloud[cloud])))
        if cloud:
            # all clouds except the base are displaced by the 
            # current translation
            print("\ttranslating cloud '%s' by " % l_cloudFile[cloud], end="")
            print(reltran)
            lC_cloud[cloud].cloud(lM_cloudOrig[cloud] + reltran)
            lv_norm[cloud] = cloud_normalize(lC_cloud[cloud].cloud())
        # We only need to process the "base" cloud ONCE in repeated
        # reltrans lookups...
        if (not b_baseProcessed and not cloud) or (b_baseProcessed and cloud):
            print("\tnormalizing...")
            lC_cloud[cloud].normalize(Gb_normalize)
            lC_cloud[cloud].usePercentiles(Gb_usePercentiles)
            lC_cloud[cloud].percentile(G_f_percentile)
            lC_cloud[cloud].asymmetricalDeviations(Gb_asymmetricalDeviations)
            lC_cloud[cloud].centerMean(Gstr_centerMean)
            print("\tfinding confidence boundary...")
            lC_cloud[cloud].confidenceBoundary_find()
            b_baseProcessed = True
        print("\textracting polygonPoints...")
        ll_polygonPoints[cloud] = lC_cloud[cloud].boundary()
        if Gb_convexHullUse:
            print("\tshaping convex hull...")
            ll_polygonPoints[cloud] = \
                convexHull_boundaryFind(ll_polygonPoints[cloud])
        print("\t...")
        lM_cloud[cloud]         = lC_cloud[cloud].cloud()
        p.append(plot(lM_cloud[cloud][:,0], lM_cloud[cloud][:,1], 
                        color=l_cloudColor[cloud], marker='*', ls='None',
                        zorder = 1))
        poly.append(deviation_plot(ll_polygonPoints[cloud], l_kernelColor[cloud])) 
        if Gb_extentReport: print(C_cloud[cloud].projectionExtent_report())

    if len(l_cloudFile) > 1:
        print("Processing statistics...")
        l_combinations = list(itertools.combinations(range(len(l_cloudFile)), 2))
        # print(l_combinations)
        for combination in l_combinations:
            g1  = combination[0]
            g2  = combination[1]
            _str_title = '%s-%s-x%5.3f-y%5.3f' % (
            		os.path.splitext(os.path.basename(l_cloudFile[g1]))[0],
            		os.path.splitext(os.path.basename(l_cloudFile[g2]))[0],
            		reltran[0], reltran[1]
            		)
            title(_str_title)
            # print("group 1 = %d" % g1)
            # print("group 2 = %d" % g2)
            v_tstat, v_pval = stats.ttest_ind(lM_cloud[g1], lM_cloud[g2], equal_var = False)
            f_tstatvNorm, f_pvalvNorm = stats.ttest_ind(lv_norm[g1], lv_norm[g2], equal_var = False)
            f_pvalMin                           = np.amin(v_pval)
            f_pvalmNorm                         = np.linalg.norm(v_pval)
            if not G_depth:
                v_ind                           = np.array( (0, 0) )
            else:
                v_ind                           = reltran/v_reltranScale + \
                                                        v_indoffset - v_nonBaseOffset
            
            M_pval[v_ind[0],v_ind[1]]           = f_pvalMin
            M_pvalmNorm[v_ind[0],v_ind[1]]      = f_pvalmNorm
            M_pvalvNorm[v_ind[0],v_ind[1]]      = f_pvalvNorm
            print("\tmin  p-val for comparison between group '%s' and '%s' is    %f" % \
                (l_cloudFile[g1], l_cloudFile[g2], f_pvalMin))
            print("\tnorm p-val for comparison between group '%s' and '%s' is    %f" % \
                (l_cloudFile[g1], l_cloudFile[g2], f_pvalmNorm))
            print("\tp-val for comparison between normed groups '%s' and '%s' is %f" % \
                (l_cloudFile[g1], l_cloudFile[g2], f_pvalvNorm))
            print("\tsaving figure '%s'..." % _str_title)
            misc.mkdir(_str_title)
            #aspectRatio_square(aspect=1)
            savefig('%s/%s.pdf' % (_str_title, _str_title), bbox_inches=0)
            savefig('%s/%s.png' % (_str_title, _str_title), bbox_inches=0)
            # Save pval matrices on each loop... allows for some data storage
            # even while processing is incomplete.
            np.savetxt('%s/%s-pval.txt'      % (_str_title, _str_title), M_pval,          fmt='%10.7f')
            np.savetxt('%s/%s-pvalmNorm.txt' % (_str_title, _str_title), M_pvalmNorm,     fmt='%10.7f')
            np.savetxt('%s/%s-pvalvNorm.txt' % (_str_title, _str_title), M_pvalvNorm,     fmt='%10.7f')
    else:
            _str_title = '%s' % (
                        os.path.splitext(os.path.basename(l_cloudFile[0]))[0],
                        )
            title(_str_title)
            misc.mkdir(_str_title)
            print("\tsaving figure '%s'..." % _str_title)
            savefig('%s/%s.pdf' % (_str_title, _str_title), bbox_inches=0)
            savefig('%s/%s.png' % (_str_title, _str_title), bbox_inches=0)

    indexTotal += 1


    # show()

