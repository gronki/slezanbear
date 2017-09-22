#!/usr/bin/env python
# coding: utf-8

#-------------------------------------------------------------------------------

from numpy.random import uniform, normal
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#-------------------------------------------------------------------------------

nx = 20
noise = 0.02
niter = 60
kmin, kmax = 0.3, 2.0

#-------------------------------------------------------------------------------

from numpy import exp
f = lambda x,k: exp(-x*k)
fk = lambda x,k: -x * exp(-x*k)
fkk = lambda x,k: x**2 * exp(-x*k)

#-------------------------------------------------------------------------------

K0 = 1.0 # uniform(kmin,kmax)
x = uniform(0,10,nx)
D = f(x,K0) + normal(0, noise, nx)

#-------------------------------------------------------------------------------

from numpy import ndarray
k = ndarray(niter)
k[0] = 2.0 # uniform(kmin,kmax)
a = 1.0
print "{:6} {:11} {:11} {:6} {:11}".format("K","Z","a","Knext","K-K0")
for i in range(1,niter):
    dY = f(x, k[i-1]) - D
    Z = sum(fk(x,k[i-1]) * dY)
    # Zk = sum(fkk(x,k[i-1]) * dY + fk(x,k[i-1])**2)
    # a = float(i) / niter
    k[i] = k[i-1] - Z * a
    print "{:6.3f} {:11.3e} {:11.3e} {:6.3f} {:11.3e}".format(k[i-1], Z, a, k[i], k[i] - K0)
    # a = (k[i] - k[i-1]) / abs(sum(fk(x,k[i]) * dY) - Z)

print "real K = {:.3f}".format(K0)

#-------------------------------------------------------------------------------

from numpy import linspace
k0 = linspace(kmin, kmax + 0.5 * (kmax - kmin), 256)
E = ndarray(k0.shape[0])
Z = ndarray(k0.shape[0])
Zk = ndarray(k0.shape[0])
for i in range(k0.shape[0]):
    dY = f(x, k0[i]) - D
    E[i] = sum(dY**2) / 2
    Z[i] = sum(fk(x,k0[i]) * dY)
    Zk[i] = sum(fkk(x,k0[i]) * dY + fk(x,k0[i])**2)

plt.figure()

ax = plt.subplot(211)
ax.set_xlim(k0.min(),k0.max())
ax.plot(k0,E)
ax.set_yscale('log')

ax = plt.subplot(212)
ax.set_xlim(k0.min(),k0.max())
ax.plot(k0,Z,'-')
ax.plot(k0,Zk,'--')
# plt.plot(k0,-Z/Zk)
ax.set_yscale('symlog', linthreshy = 0.1)
for i in range(niter):
    ax.axvline(k[i], color = cm.inferno(1-float(i) / niter))
ax.axvline(K0, linewidth = 1.6, linestyle = '--', color = 'black')
ax.grid()

plt.savefig('/tmp/optim2k.png')
# plt.show()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


x0 = linspace(0,x.max(),200)
plt.figure()
plt.scatter(x,D)
plt.plot(x0, f(x0,K0), color = '#9C9C9C')
for i in range(niter): plt.plot(x0, f(x0,k[i]), color = cm.inferno(1-float(i) / niter))
plt.savefig('/tmp/optim2.png')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
