#!/usr/bin/env python
#
# Plotting program, to use to visualize the output of the mini-app
# CSCS - Summer School 2015
# Written by Jean M. Favre
#
# updated to Python3 (Feb 2020)

# check the DISPLAY environment variable
# if it is not defined, don't show the image
import os
if os.environ.get('DISPLAY') == None:
    import matplotlib
    matplotlib.use('Agg')

import pylab as pl
import numpy as np

# change filename below as appropriate
headerFilename =  "output.bov"
headerfile = open(headerFilename, "r")
header = headerfile.readlines()
headerfile.close()

rawFilename = header[1].split()[1]

res_x, res_y = [int(x) for x in header[2].strip().split()[1:3]]

# for python 2.7
#print 'xdim %s' % res_x
#print 'ydim %s' % res_y

# for python3
print('xdim %s' % res_x)
print('ydim %s' % res_y)

data = np.fromfile(rawFilename, dtype=np.double, count=res_x*res_y, sep='')
assert data.shape[0] == res_x * res_y, "raw data array does not match the resolution in the header"

x = np.linspace(0., 1., res_x)
y = np.linspace(0., 1.*res_y/res_x, res_y)
X, Y = np.meshgrid(x, y)
# number of contours we wish to see
N = 12
pl.contourf(X, Y, data.reshape(res_y, res_x), N, alpha=.75, cmap='jet')
C=pl.contour(X, Y, data.reshape(res_y, res_x), N, colors='black', linewidths=.1)

pl.clabel(C, inline=1)
pl.axes().set_aspect('equal')
pl.savefig("output.png", dpi=72)
pl.show()
