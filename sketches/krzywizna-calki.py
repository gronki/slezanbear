#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import latexrc
from numpy import linspace
from numpy import sqrt, cos, arccos, exp

#-------------------------------------------------------------------------------

# parameters
R = 6371.0
H = R / 50
Hscale = R / 500
Dmax = R * arccos(R / (R + H))
D = Dmax * 0.7
n = 4000

print "R = {:.2f}".format(R)
print "H = {:.2f}".format(H)
print "D = {:.2f}".format(D)
print "Hscale = {:.2f}".format(Hscale)

#-------------------------------------------------------------------------------

clr_flt = '#A5AF3A'
clr_rnd = '#0C62D5'
clr_apx = '#FF1F00'

from numpy import diff,cumsum
midpoints = lambda x: ( x[1:] + x[:-1] ) / 2
integral = lambda y,x: cumsum( midpoints(y) * diff(x) )

#-------------------------------------------------------------------------------

print u'FLAT EARTH'

L0 = sqrt(D**2 + H**2)
l0 = linspace(0,L0,n)
h0 = (H / L0) * l0

print "L = {:.2f}".format(L0)

ch0 = exp(-h0 / Hscale)
tau0 = integral(ch0,l0)
tau0th = Hscale * (1 - exp(- (H / L0) * l0 / Hscale)) / (H / L0)

#-------------------------------------------------------------------------------

fig,ax = plt.subplots()
ax.set_xlabel('$l$')
ax.set_ylabel('$\\tau$')
ax.set_title('Flat Earth')
ax.plot(l0[1:], tau0, color = '#19B2E3', label = 'numerical')
ax.plot(l0, tau0th, '--', color = '#103885', label = 'analytical')
ax.legend(loc = 'lower right')

#-------------------------------------------------------------------------------

print u'ROUND EARTH'

L1 = sqrt(H**2 + 2*R*(R+H)*(1 - cos(D/R)))
l1 = linspace(0,L1,n)
A = H * (H + 2*R) / L1**2
h1 = sqrt(R**2 + l1*L1*(A-1) + l1**2) - R

print "L = {:.2f}".format(L1)
print "A = {:.2f}".format(A)

ch1 = exp(-h1 / Hscale)
tau1 = integral(ch1,l1)

#-------------------------------------------------------------------------------

print u'APPROXIMATION'

L2 = sqrt(H**2 + (1 + H/R) * D**2)
l2 = linspace(0,L2,n)
h2 = l2 * (H / L2 - (L2 - l2) / (2*R))

print "L = {:.2f}".format(L2)

ch2 = exp(-h2 / Hscale)
tau2 = integral(ch2,l2)

from numpy import pi
from erf import erf

def liczcalke(R,D,H,L,l):
    A = 2*R*H / L**2
    D = sqrt(2*R*Hscale) / L
    F = (A - 1) / (2*D)
    f = l / sqrt(2*R*Hscale)
    return L * sqrt(pi/4) * exp(F**2) * D * ( erf(F + f) - erf(F) )

tau2th = liczcalke(R,D,H,L2,l2)

#-------------------------------------------------------------------------------

fig,ax = plt.subplots()
ax.set_xlabel('$l$')
ax.set_ylabel('$\\tau$')
ax.set_title('Approximation Earth')
ax.plot(l2[1:], tau2, color = '#19B2E3', label = 'numerical')
ax.plot(l2, tau2th, '--', color = '#103885', label = 'analytical')
ax.legend(loc = 'lower right')

plt.savefig('/tmp/przyblizenie.pdf')

#-------------------------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize = (5.9, 3.6))

ax = axes[0]
ax.set_xlabel('$l$')
ax.set_ylabel('$h(l)$')
ax.set_title('height vs path -- all three')
ax.plot(l0, h0, color = clr_flt, label = 'Flat Earth')
ax.plot(l1, h1, color = clr_rnd, label = 'Round Earth')
ax.plot(l2, h2, color = clr_apx, label = 'Approximation')
ax.legend(loc = 'lower right')

ax = axes[1]
ax.set_xlabel('$l$')
ax.set_ylabel('$\\tau$')
ax.set_title('flat vs round earth')
ax.plot(l0[1:], tau0, color = clr_flt, label = 'Flat Earth')
ax.plot(l1[1:], tau1, color = clr_rnd, label = 'Round Earth')
ax.plot(l2[1:], tau2, color = clr_apx, label = 'Approximation')
ax.legend(loc = 'lower right')

plt.subplots_adjust(0.08,0.12,0.95,0.92,0.25,0.25)

plt.savefig('/tmp/krzywizna.pdf')

#-------------------------------------------------------------------------------

from numpy import array
d = linspace(0,Dmax,n)

L0 = sqrt(H**2 + d**2)
tau0 = L0 * (Hscale / H) * (1 - exp(-H / Hscale))

L1 = sqrt(H**2 + (1 + H/R) * d**2)
A = 2*R*H / L1**2
B = sqrt(2*R*Hscale) / L1
tau1 = L1 * sqrt(pi/4) * exp((A - 1)**2 / (2*B)**2) * B \
    * ( erf((A + 1) / (2*B)) - erf((A - 1) / (2*B)) )
tau2 = L1 * exp((1 - 2*A) / (2*B)**2)
tau2 = L1 * exp(L1**2 / (8 * R * Hscale)) * exp(- 0.5 * H / Hscale)

kkk = exp(- (A + 0.8)**2 / (2*B)**2) \
    + exp(- (A + 0.4)**2 / (2*B)**2) \
    + exp(- (A + 0.0)**2 / (2*B)**2) \
    + exp(- (A - 0.4)**2 / (2*B)**2) \
    + exp(- (A - 0.8)**2 / (2*B)**2)
tau3 = L1 * (kkk / 5) * exp((A - 1)**2 / (2*B)**2)

aa = lambda l: exp(-l / Hscale * ( H/L1 - (L1-l)/(2*R)))
tau4 = aa(0.1*L1) + aa(0.3*L1) + aa(0.5*L1) + aa(0.7*L1) + aa(0.9*L1)
tau4 *= L1 / 5

fig,axes = plt.subplots(2,2,figsize = (9,9))
ax = axes[0,0]
ax.plot(d,tau0, color = '#D7B212', label = 'flat')
ax.plot(d,tau1, color = '#006AB6', label = 'exact')
ax.plot(d,tau2, '--', color = '#31BD00', label = 'appx')
ax.plot(d,tau3, '--', color = '#CF03DA', label = 'appx2')
ax.plot(d,tau4, linewidth = 0.5, color = '#5E0091', label = 'appx3')
ax.set_ylabel('$\\tau')
# ax.set_ylim(0,2 * Hscale / H)
ax.legend(loc = 'best', fontsize = 10)

ax = axes[1,0]
# ax.plot(d, (A-1) / (2*B))
# ax.plot(d, A / (2*B))
# ax.plot(d, (A+1) / (2*B))
# ax.plot(d, 2 / (2*B), ':')
ax.plot(d,A)

ax = axes[0,1]

ax.plot(d, exp(-(A - 1)**2 / (2*B)**2), ':', label = '$\\exp\\left(-\\frac{(A - 1)^2}{4B^2}\\right)$')
ax.plot(d, sqrt(pi/4)*B*( erf((A + 1) / (2*B)) - erf((A - 1) / (2*B)) ), label = '$\\sqrt{\\pi/4}B\\left( {\\rm erf}\\left(\\frac{A + 1}{2B}\\right) - {\\rm erf}\\left(\\frac{A - 1}{2B}\\right) \\right)$')
ax.plot(d, exp(- A**2 / (2*B)**2), '--', label = '$\\exp\\left(- \\frac{A^2}{4B^2}\\right)$')
ax.plot(d, exp(- (2*A+1)**2 / (4*B)**2) + exp(- (2*A-1)**2 / (4*B)**2), color = '#CF03DA', linestyle = ':', label = 'new')
ax.set_yscale('log')
ax.legend(loc = 'best', fontsize = 10)

plt.savefig('/tmp/krzywizna2.pdf')

#-------------------------------------------------------------------------------

# plt.show()
