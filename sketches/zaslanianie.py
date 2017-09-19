from numpy import linspace, logspace, meshgrid, where, arange, interp, sign
from numpy import sin, cos, log, sqrt, arccos, exp, abs, pi, arcsinh, log10
from numpy.random import uniform
import matplotlib.pyplot as plt
# import latexrc

#-------------------------------------------------------------------------------

xmax = 30.0
hmax = 10.0
x = linspace(-xmax, xmax, 400)
dx = abs(x[1] - x[0])
h0 = linspace(0, hmax, 200)

D,Habs = meshgrid(x, h0, indexing = 'xy')

#-------------------------------------------------------------------------------

Hterr = x * 0
for i in range(int(xmax * 3)):
    w = 10**uniform(log10(0.1),log10(4.0))
    xc = uniform(-(xmax + 2*w), xmax + 2*w)
    h = 10**uniform(-2.0, 0.0)
    Hterr += h * sqrt(w) * exp(- 0.5 * (x - xc)**2 / w**2)

Hsrc = interp(0, x, Hterr)
H = Habs - Hsrc

#-------------------------------------------------------------------------------

def diff2(x,y):
    from numpy import ndarray
    assert x.shape == y.shape
    d = ndarray(y.shape)
    d[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    d[0] = (y[1] - y[0]) / (x[1] - x[0])
    d[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    return d

tandel = diff2(x, Hterr)
cosdel = tandel / sqrt(1 + tandel**2)

#-------------------------------------------------------------------------------

Re = 6341.0
R = Re + Hsrc
L = sqrt(H**2 + (1 + H / R) * D**2)
Q = D**2 / (2*R*H)
costh = H / L * (1 + Q)
cosgam = H / L * (1 - Q * (1 + H / R))

I = 1 / ((2 * L / xmax)**2 + 1)
from threshold import th

for j in range(H.shape[1]):
    for i in range(H.shape[0]):
        lray = linspace(0, L[i,j], max(L[i,j] / dx,10.0))
        hray = lray * (H[i,j] / L[i,j] - (L[i,j] - lray) / (2 * R**2))
        dray = sqrt(lray**2 - hray**2) / sqrt(1 + hray / R)
        hter = interp(sign(x[j]) * dray, x, Hterr - Hsrc)
        I[i,j] *= (2 + th(Habs[i,j] - Hterr[j], 0.03)) / 3
        I[i,j] *= (2 + th((hray - hter).min(), 0.03)) / 3

#-------------------------------------------------------------------------------

fig,ax = plt.subplots( figsize = (18,5))

ax.pcolor(x, h0, I, cmap = 'Blues_r')
ax.plot(x, Hterr, color = 'black', label = 'terrain')
ax.set_xlim(-xmax,xmax)
ax.set_aspect(1)

plt.tight_layout()

plt.savefig('/tmp/zaslanianie.png')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
