from numpy import linspace, meshgrid, ndarray, where, sum, pi, abs
from numpy import cos, sin, exp, log, log10, sqrt, arcsinh, arctan

r = 1
y = linspace(0,r,120)
x = linspace(-r,r,2*len(y))

xx, yy = meshgrid(x,y)

rr = sqrt(xx**2 + yy**2)
costh0 = yy / rr

A0 = pi

a = linspace(-0.5,0.5)
costh1 = ndarray(xx.shape)

def fun(x,z):
    plx = linspace(-0.5,0.5,51)
    ply = linspace(-0.5,0.5,51)
    plxx, plyy = meshgrid(plx, ply)
    rr = sqrt(plxx**2 + plyy**2)
    mask = where(rr <= 0.5, 1, 0) * (1 - (rr / 0.5)**2)
    plxx += x
    rrr = sqrt(plxx**2 + plyy**2 + z**2)
    z2 = z / rrr
    da = abs((plx[1] - plx[0])*(ply[1] - ply[0]))
    z1 = da / (da + 2 * pi * rrr**2) * (z / rrr)
    return sum(z2 * z1 * mask) / sum(z1 * mask)

for i in range(xx.shape[0]):
    for j in range(xx.shape[1]):
        costh1[i,j] = fun(xx[i,j], yy[i,j])

# ztan = lambda x: x/sqrt(x**2 + 1)
costh2 = yy * ( arcsinh((xx + 0.5) / yy) - arcsinh((xx - 0.5) / yy) )
# costh2 =  ( arctan ((xx + 0.5) / yy) - arctan ((xx - 0.5) / yy) ) \
#         / ( arcsinh((xx + 0.5) / yy) - arcsinh((xx - 0.5) / yy) )

# aa = exp( -0.5 * rr**2 / 0.5**2)
# aa = where(rr > 0.5, 0, 1)
# a = A0 / (A0 + 2*pi*rr**2)
# costh2 = yy / rr * a + 0.5 * (1-a)

import matplotlib.pyplot as plt

fig, axes = plt.subplots(2,2,figsize=(16,9))

c_levels = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 1.0]

ax = axes[0,0]
m = ax.contourf(x, y, costh0, levels = c_levels, cmap = 'inferno')
plt.colorbar(m, ax = ax)
ax.set_aspect('equal')

ax = axes[1,0]
m = ax.contourf(x, y, costh1, levels = c_levels, cmap = 'inferno')
plt.colorbar(m, ax = ax)
ax.set_aspect('equal')

ax = axes[1,1]
m = ax.contourf(x, y, costh2, levels = c_levels, cmap = 'inferno')
plt.colorbar(m, ax = ax)
ax.set_aspect('equal')

ax = axes[0,1]
m = ax.pcolor(x, y, costh2 / costh1 - 1, vmin = -0.5, vmax = 0.5, cmap = 'seismic')
ax.set_aspect('equal')
plt.colorbar(m, ax = ax)

plt.show()
