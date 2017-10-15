import matplotlib.pyplot as plt
from numpy import pi, arccos, cos, sin, sqrt, log, log10, exp
from numpy import linspace, meshgrid, where, ndarray, array, float32, float64
from slezanbear import sbmap
from threshold import th

radearth = 6371 * 1e3

dist = lambda f1,l1,f2,l2: arccos( \
        cos(f1 * (pi / 180)) * cos(f2 * (pi / 180)) \
    *   cos((l1 - l2) * (pi / 180)) \
    +   sin(f1 * (pi / 180)) * sin(f2 * (pi / 180)))

latlim  = 50.75, 51.05
lnglim  = 16.55, 16.85

nord = 5

mkgt = lambda lat,lng: array([lng[0], lng[1] - lng[0], 0,
            lat[0], 0, lat[1] - lat[0]],
            dtype = float32, order = 'F')

#-------------------------------------------------------------------------------

def latlng(latlim, lnglim, nx, ny):
    lat0 = linspace(latlim[0], latlim[1], ny)
    lng0 = linspace(lnglim[0], lnglim[1], nx)
    lat, lng = meshgrid(lat0, lng0, indexing = 'ij')
    gt = array([ lnglim[0], (lnglim[1] - lnglim[0]) / (nx - 1), 0,
                 latlim[0], 0, (latlim[1] - latlim[0]) / (ny - 1)],
                 dtype = float32, order = 'F')
    return lat, lng, gt

#-------------------------------------------------------------------------------

gauss = lambda x,w: exp(-0.5 * x**2 / w**2)

n = 2**(nord + 1)
lat, lng, gth = latlng(latlim, lnglim, n, n)
r1 = radearth * dist(50.860, 16.740, lat, lng)
r2 = radearth * dist(50.875, 16.703, lat, lng)
mountain = lambda f,l,h,w: h * exp(- 0.5 * radearth**2 * dist(f,l,lat,lng)**2 / w**2)
hterr = array(200 + 300 * th(3000 - r1, 1200) + 400 * th(1800 - r2, 800) \
        - 300 * (lat - 50.90) - 30 * (lng - 16.70) \
        + mountain(50.951, 16.700, 500, 500) \
        + mountain(50.951, 16.710, 500, 200) \
        + mountain(50.955, 16.720, 400, 400) \
        + mountain(50.955, 16.725, 300, 600) \
        + mountain(50.960, 16.730, 400, 500) \
        + mountain(50.891, 16.655, 500, 1300),
    copy = True, dtype = float32, order = 'F')

#-------------------------------------------------------------------------------

n = 2**(nord + 1)
lat, lng, gti = latlng(latlim, lnglim, n, n)
r = radearth * dist(50.917, 16.720, lat, lng)
I0 = array(where(r < 1e3 , 99, 0),
    copy = True, dtype = float32, order = 'F')

#-------------------------------------------------------------------------------

n = 2**nord
lat, lng, gt = latlng(latlim, lnglim, n, n)
I1 = ndarray((n,n), dtype = float64, order = 'F')
I2 = ndarray((n,n), dtype = float64, order = 'F')

#-------------------------------------------------------------------------------


sbmap(  I0,     I0.shape[0],    I0.shape[1],    gti,
        hterr,  hterr.shape[0], hterr.shape[1], gth,
        I1, I2, I1.shape[0],    I1.shape[1],    gt)

fig, axes = plt.subplots(2, 2, figsize = (13,11))
p = axes[0,0].imshow(hterr, interpolation = 'none', cmap = 'BrBG_r',
    extent = [lnglim[0], lnglim[1], latlim[0], latlim[1]], origin = 'lower')
plt.colorbar(p, ax = axes[0,0])

p = axes[0,1].imshow(I0, interpolation = 'none', cmap = 'hot',
    extent = [lnglim[0], lnglim[1], latlim[0], latlim[1]], origin = 'lower')
plt.colorbar(p, ax = axes[0,1])

p = axes[1,0].imshow(sqrt(I2), interpolation = 'none', cmap = 'inferno',
    extent = [lnglim[0], lnglim[1], latlim[0], latlim[1]], origin = 'lower')
plt.colorbar(p, ax = axes[1,0])

p = axes[1,1].imshow(I2 / I1 - 1, interpolation = 'none', cmap = 'inferno',
    extent = [lnglim[0], lnglim[1], latlim[0], latlim[1]], origin = 'lower', vmax = 0)
plt.colorbar(p, ax = axes[1,1])

print '/tmp/testcase1.png'
plt.savefig('/tmp/testcase1.png')
# plt.show()
