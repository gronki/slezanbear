# coding: utf-8
import gdal
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from re import match, sub

i2m = lambda i: 5.5 - 2.5 * np.log10(i)
m2i = lambda m: 10 ** ((5.5 - m) / 2.5)

colors = [
    (18.000, '#FFFFFF'),
    (19.400, '#D2613D'),
    (20.000, '#FFDE00'),
    (20.800, '#476100'),
    (21.100, '#2E98A6'),
    (21.450, '#0E0036'),
    (21.600, '#3A0047'),
    (21.700, '#1C1C1C'),
    (22.000, '#000000'),
]

# ccolors = [ c for l,c in colors ]
# clevels = [ l for l,c in colors ]

sg = { 'red': [], 'green': [], 'blue': [] }
for l,c in colors:
    cr = int(c[1:3], 16)
    cg = int(c[3:5], 16)
    cb = int(c[5:7], 16)
    sg['red'  ].append(( (l - 18.0) / (22.0 - 18.0), cr / 255.0, cr / 255.0))
    sg['green'].append(( (l - 18.0) / (22.0 - 18.0), cg / 255.0, cg / 255.0))
    sg['blue' ].append(( (l - 18.0) / (22.0 - 18.0), cb / 255.0, cb / 255.0))

clevels = np.concatenate([
    np.arange(18.0, 20.0, 0.250),
    np.arange(20.0, 20.8, 0.125),
    np.arange(20.8, 21.5, 0.050),
    np.arange(21.5, 21.8, 0.025),
    np.array([22.0]) ])

from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap('sky', sg)

for fn in argv[1:]:

    if not match(r'.*\.tif$', fn): continue

    gd = gdal.Open(fn)
    gt = gd.GetGeoTransform()
    nx = gd.RasterXSize
    ny = gd.RasterYSize
    arr = gd.ReadAsArray()
    hobs = arr[0,:,:]
    sky = arr[1:,:,:]

    fs = ( 10 * np.sqrt(nx/float(ny)), 10 / np.sqrt(nx/float(ny)) )

    ex = [ gt[0], gt[0] + nx * gt[1],
           gt[3] + ny * gt[5], gt[3] ]


    x = np.linspace(gt[0], gt[0] + nx * gt[1], nx)
    y = np.linspace(gt[3], gt[3] + ny * gt[5], ny)

    #---------------------------------------------------------------------------

    fig,ax = plt.subplots(figsize = fs)
    m = ax.imshow(hobs, cmap = 'BrBG_r')
    ax.set_aspect('auto')

    plt.colorbar(m, ax = ax)
    plt.savefig(sub(r'\.tif$','.el.png',fn), dpi = 144, interpolation = 'none')

    #---------------------------------------------------------------------------

    fig,ax = plt.subplots(figsize = fs)
    m = ax.contourf(x, y, i2m(sky[3,:,:]), levels = clevels, cmap = cmap)
    plt.colorbar(m, ax = ax)
    ax.set_aspect('auto')

    plt.savefig(sub(r'\.tif$','.png',fn), dpi = 144)

    #---------------------------------------------------------------------------

    fig,ax = plt.subplots(figsize = fs)
    m = ax.contourf(x, y, i2m(sky[2,:,:]), levels = clevels, cmap = cmap)
    # m = ax.imshow(i2m(sky[2,:,:]), cmap = cmap, vmin = 18.0, vmax = 22.0)
    ax.set_aspect(1 / np.cos(51.0 * np.pi / 180))
    ax.set_aspect('auto')
    plt.colorbar(m, ax = ax)
    plt.savefig(sub(r'\.tif$','.noat.png',fn), dpi = 144)

    #---------------------------------------------------------------------------

    fig,ax = plt.subplots(figsize = fs)
    m = ax.imshow(-2.5 * np.log10(sky[3,:,:] / sky[2,:,:]), cmap = 'inferno_r',  vmin = 0.0, vmax = 0.1, interpolation = 'none')
    ax.set_aspect('auto')
    plt.colorbar(m, ax = ax)
    plt.savefig(sub(r'\.tif$','.atdiff.png',fn), dpi = 144)

    #---------------------------------------------------------------------------

    fig,ax = plt.subplots(figsize = fs)
    m = ax.imshow(-2.5 * np.log10(sky[2,:,:] / sky[1,:,:]), cmap = 'RdGy',  vmin = -0.1, vmax = 0.1, interpolation = 'none')
    ax.set_aspect('auto')
    plt.colorbar(m, ax = ax)
    plt.savefig(sub(r'\.tif$','.phdiff.png',fn), dpi = 144)
