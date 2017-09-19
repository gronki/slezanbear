#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import latexrc
from numpy import linspace
from numpy import sqrt, cos, arccos, exp

#-------------------------------------------------------------------------------

# parameters
R = 6371.0
H = R / 10
Dfrac = 0.999
n = 4000

#-------------------------------------------------------------------------------

def main(R, H, Dfrac, ax):

    #---------------------------------------------------------------------------

    Dmax = R * arccos(R / (R + H))
    D = Dmax * Dfrac

    #---------------------------------------------------------------------------

    print "R = {:.2f}".format(R)
    print "H = {:.2f}".format(H)
    print "D = {:.2f}".format(D)

    #---------------------------------------------------------------------------

    clr_flt = '#A5AF3A'
    clr_rnd = '#0C62D5'
    clr_apx = '#FF1F00'

    #---------------------------------------------------------------------------

    print u'FLAT EARTH'

    L0 = sqrt(D**2 + H**2)
    l0 = linspace(0,L0,n)
    h0 = (H / L0) * l0

    print "L = {:.2f}".format(L0)

    #---------------------------------------------------------------------------

    print u'ROUND EARTH'

    L1 = sqrt(H**2 + 2*R*(R+H)*(1 - cos(D/R)))
    l1 = linspace(0,L1,n)
    A = H * (H + 2*R) / L1**2
    h1 = sqrt(R**2 + l1*L1*(A-1) + l1**2) - R

    print "L = {:.2f}".format(L1)
    print "A = {:.2f}".format(A)


    #---------------------------------------------------------------------------

    print u'APPROXIMATION'

    L2 = sqrt(H**2 + (1 + H/R) * D**2)
    l2 = linspace(0,L2,n)
    h2 = l2 * (H / L2 - (L2 - l2) / (2*R))

    print "L = {:.2f}".format(L2)

    #---------------------------------------------------------------------------

    print u'APPROXIMATION'

    L3 = sqrt(H**2 + (1 + H/R) * D**2)
    l3 = linspace(0,L3,n)
    h3 = l3 * (H / L3) * (1 - (L3**2 - H**2 - L3*l3) / (2*R*H))

    print "L = {:.2f}".format(L3)

    #---------------------------------------------------------------------------


    clr = ['#F5B920', '#003173', '#51C709', '#FF2DCA']

    ax.set_title('$H/R = 1/{}$, $D / D_{{\\rm horiz}} = {}$'.format(int(round(R/H)), Dfrac))

    ax.grid()

    ax.plot(l0 / L0, h0 / H, color = clr[0], label = 'flat earth')
    ax.plot(l2 / L0, h2 / H, color = clr[2], label = 'approximation')
    ax.plot(l3 / L0, h3 / H, color = clr[3], linewidth = 0.7, label = 'approximation 2')
    ax.plot(l1 / L0, h1 / H, color = clr[1], linestyle = '--', linewidth = 1.2, label = 'round earth (exact)')

    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel('$l / \\sqrt{H^2+D^2}$')
    ax.set_ylabel('$h(l) / H$')

#-------------------------------------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize = (6,6))
main(R, R / 2, 0.25, axes[0,0])
main(R, R / 2, 0.99, axes[0,1])
main(R, R / 20, 0.25, axes[1,0])
main(R, R / 20, 0.99, axes[1,1])
axes[1,0].legend()
plt.tight_layout()
plt.savefig('/tmp/wysokosc.pdf')
