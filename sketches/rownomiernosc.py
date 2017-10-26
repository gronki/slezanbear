
from numpy import linspace, meshgrid, arange
from numpy import sin, sqrt, log, log10, exp, abs, pi

rnge = 15
nodx, nody = meshgrid(arange(-rnge,rnge) + 0.5, arange(-rnge,rnge) + 0.5)
A = 1.0

x, y = meshgrid(linspace(-rnge,rnge,300), linspace(-rnge,rnge,300))

sky = x * 0
H = sin( pi / 2 * x / rnge ) + 1

#--------------------------------------------------------------------------

for i in range(nodx.shape[0]):
    for j in range(nodx.shape[1]):
        D = sqrt((x - nodx[i,j])**2 + (y - nody[i,j])**2)
        L = sqrt((x - nodx[i,j])**2 + (y - nody[i,j])**2 + H**2)
        cosg = H / L
        K = pi * L*H / A
        sky += 2 * pi * cosg * A / (2 * pi * L**2 + cosg * A)


#--------------------------------------------------------------------------

import matplotlib.pyplot as plt
fig,ax = plt.subplots(figsize = (9,9))
ax.imshow(sky / (2*pi) - 1, extent = [-rnge, rnge, -rnge, rnge], vmin = -1, vmax = 1, cmap = 'seismic', interpolation = 'bilinear')
plt.savefig('/tmp/c.png')
