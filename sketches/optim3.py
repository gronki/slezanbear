#!/usr/bin/env python
# coding: utf-8

#-------------------------------------------------------------------------------

from numpy import exp, sqrt, linspace, logspace, pi, log
from numpy.random import uniform, normal
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#-------------------------------------------------------------------------------

G = lambda x,a,b,c: b + (a-b)*exp(-x/c)
Ga = lambda x,a,b,c: exp(-x/c)
Gb = lambda x,a,b,c: 1 - exp(-x/c)
Gc = lambda x,a,b,c: x*(a - b)*exp(-x/c)/c**2

#-------------------------------------------------------------------------------

npoints = 70
x0 = uniform(0,15,npoints)
a0 = uniform(0,10)
b0 = a0 + uniform(-10,10)
c0 = uniform(2,8)
y0 = G(x0,a0,b0,c0) + normal(0, 10**uniform(-3,-0.5), npoints)

print a0,b0,c0

#-------------------------------------------------------------------------------

from numpy import zeros, array
niter = 48
P = zeros((3,niter))
P[:,0] = [0,0,5]
for i in range(1,niter):
    F = 0.5 * sum((G(x0, *P[:,i-1]) - y0)**2)
    dy = G(x0, *P[:,i-1]) - y0
    da = sum(Ga(x0, *P[:,i-1]) * dy)
    db = sum(Gb(x0, *P[:,i-1]) * dy)
    dc = sum(Gc(x0, *P[:,i-1]) * dy)
    alp = 1.0
    for j in range(30):
        Ptest = P[:,i-1] - array([da,db,dc]) * alp
        Ftest = 0.5 * sum((G(x0, *Ptest) - y0)**2)
        print 'alp = {:8.3f}    Ftest = {:12.4e}'.format(alp, Ftest/F - 1)
        if Ftest < F: break
        alp *= 0.5
    P[:,i] = Ptest
    print " {:8.3f} {:8.3f} {:8.3f}".format(*P[:,i])
    print "({:8.3f} {:8.3f} {:8.3f})".format(a0,b0,c0)

#-------------------------------------------------------------------------------

plt.figure(figsize = (12,8))
from matplotlib.cm import inferno,hot,jet,Greys
x = linspace(0,15,200)
for i in range(niter):
    plt.plot(x, G(x,*P[:,i]), color = inferno(1 - float(i) / niter))
plt.scatter(x0,y0)
plt.savefig('/tmp/optim3a.png')
#-------------------------------------------------------------------------------

plt.figure()
plt.plot(P.transpose())
plt.axhline(a0)
plt.axhline(b0)
plt.axhline(c0)
plt.savefig('/tmp/optim3.png')

# plt.scatter(x0,y0)
# x = linspace(0,15,100)
# plt.plot(x,G(x,a0,b0,c0), color = 'black', linewidth = 1.5)
# plt.show()
