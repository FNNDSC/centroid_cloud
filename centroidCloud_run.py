#!/usr/bin/env python

Gstr_matrixType = 'lower'
Gb_saveFig      = False

from    pylab import *
import  numpy as np
import  getopt
import  sys


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

M_cNr = np.genfromtxt('iter_c1_u1-Normal-r-%s.crd' % Gstr_matrixType)
M_cNg = np.genfromtxt('iter_c1_u1-Normal-g-%s.crd' % Gstr_matrixType)
M_cNb = np.genfromtxt('iter_c1_u1-Normal-b-%s.crd' % Gstr_matrixType)

M_cPr = np.genfromtxt('iter_c1_u1-PMG-r-%s.crd' % Gstr_matrixType)
M_cPg = np.genfromtxt('iter_c1_u1-PMG-g-%s.crd' % Gstr_matrixType)
M_cPb = np.genfromtxt('iter_c1_u1-PMG-b-%s.crd' % Gstr_matrixType)

figure()
#axis([0.2, 0.6, 0.4, 0.8])
grid()

p1, = plot(M_cNr[:,0], M_cNr[:,1], color='#FF0000', marker='*', ls='None');
p2, = plot(M_cNg[:,0], M_cNg[:,1], color='#00FF00', marker='*', ls='None');
p3, = plot(M_cNb[:,0], M_cNb[:,1], color='#0000FF', marker='*', ls='None');

nRr = deviation_plot(M_cNr, '#FF0000')
nRg = deviation_plot(M_cNg, '#00FF00')
nRb = deviation_plot(M_cNb, '#0000FF')

plot(M_cPr[:,0], M_cPr[:,1], color='#770000', marker='+', ls='None');
plot(M_cPg[:,0], M_cPg[:,1], color='#007700', marker='+', ls='None');
plot(M_cPb[:,0], M_cPb[:,1], color='#000077', marker='+', ls='None');

pRr = deviation_plot(M_cPr, '#770000')
pRg = deviation_plot(M_cPg, '#007700')
pRb = deviation_plot(M_cPb, '#000077')

xlabel('X-centroid position')
ylabel('Y-centroid position')
l1 = legend([nRr, nRg, nRb], ['red - normal', 'green - normal', 'blue - normal'],
        loc = 3)
l2 = legend([pRr, pRg, pRb], ['red - PMG', 'green - PMG', 'blue - PMG'], 
        loc = 4)
gca().add_artist(l1)
title('Centroid distributions for the %s matrix after evolution' % Gstr_matrixType)

if Gb_saveFig:
    print "Saving figure to png and eps...", 
    savefig('centroids-%s.png' % Gstr_matrixType)
    savefig('centroids-%s.eps' % Gstr_matrixType)
    print "Done."

show()

#raw_input('Press any key to exit...\n')


