from numpy import linspace, logspace
from numpy import sin, cos, log, sqrt, arccos, exp, abs, diff, pi
import matplotlib.pyplot as plt
import latexrc

#-------------------------------------------------------------------------------

a = linspace(0,90,9000)
x = cos(pi * a / 180)
cutoff = 3

#-------------------------------------------------------------------------------

th = lambda x: 1 / ( 1 + exp(-4*x) )
y = th(2 * (x / cos(pi * (90 - cutoff) / 180) - 1))

#-------------------------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize = (6,3))
axes[0].plot(a,y)
axes[0].set_xlabel('degrees')
axes[1].plot(x,y)
axes[1].set_xlabel('cosine')

#-------------------------------------------------------------------------------

plt.savefig('/tmp/odciecie.png')
# plt.show()

#-------------------------------------------------------------------------------
