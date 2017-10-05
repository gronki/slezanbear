#!/usr/bin/env python
# coding: utf-8

#-------------------------------------------------------------------------------

from numpy import exp, sqrt, linspace, logspace, pi, log, sin, cos, arccos, arctan2
from numpy.random import uniform, normal
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#-------------------------------------------------------------------------------

lat1,lng1 = 51.0, 17.0
lat2,lng2 = 52.0, 19.0
Re = 6371.0
Re = 180 / pi

lat1, lng1, lat2, lng2 = [ pi * x / 180 for x in [lat1, lng1, lat2, lng2]]

#-------------------------------------------------------------------------------

t = linspace(0,1,200)

#-------------------------------------------------------------------------------

Lf = Re * sqrt((lat2 - lat1)**2 + cos((lat1+lat2) / 2)**2 * (lng2 - lng1)**2)
ff = lat1 + (lat2 - lat1) * t
lf = lng1 + (lng2 - lng1) * t
print "L(flat) = {}".format(Lf)

#-------------------------------------------------------------------------------

Ls = Re * arccos(cos(lat1)*cos(lat2)*cos(lng1-lng2) + sin(lat1)*sin(lat2))
print "L(sphere) = {}".format(Ls)

from numpy import array
geo2xyz = lambda lat,lng: array([cos(lat)*cos(lng), cos(lat)*sin(lng), sin(lat)])
xyz2geo = lambda x,y,z: array([arctan2(z, sqrt(x**2 + y**2)), arctan2(y, x)])
x1 = geo2xyz(lat1,lng1)
x2 = geo2xyz(lat2,lng2)

from numpy import cross, dot, ndarray
from numpy.linalg import norm
w = cross(cross(x1,x2) / norm(cross(x1,x2)),x1)
print sum(x1**2), sum(x2**2), sum(w**2), sum(cross(x1,x2)**2)
print dot(x1,w), dot(x2,w)

fs = ndarray(t.shape[0])
ls = ndarray(t.shape[0])
for i in range(t.shape[0]):
    xi = x1 * cos(t[i] * arccos(dot(x1,x2))) + w * sin(t[i] * arccos(dot(x1,x2)))
    print 100.0 * i / t.shape[0], xi - x1, x2 - xi
    fs[i], ls[i] = xyz2geo(xi[0], xi[1], xi[2])


#-------------------------------------------------------------------------------

fig,ax = plt.subplots(figsize = (20,10))

ax.scatter([lng1* 180 / pi,lng2* 180 / pi], [lat1* 180 / pi,lat2* 180 / pi])
ax.plot(lf * 180 / pi, ff * 180 / pi, color = '#8ED6A2')
ax.plot(ls * 180 / pi, fs * 180 / pi, color = '#00B933')
# ax.set_xlim(-180,180)
# ax.set_ylim(-90,90)
ax.set_aspect(1)

#-------------------------------------------------------------------------------

plt.show()

#-------------------------------------------------------------------------------
