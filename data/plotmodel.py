# coding: utf-8
import gdal
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from re import match, sub

i2m = lambda i: 5.5 - 2.5 * np.log10(i)
m2i = lambda m: 10 ** ((5.5 - m) / 2.5)

clevels = [
    17.000,
    18.000,
    18.500,
    19.000,
    19.500,
    20.000,
    20.500,
    20.750,
    21.000,
    21.125,
    21.250,
    21.375,
    21.500,
    21.550,
    21.600,
    21.650,
    21.700,
]

ccolors = [
    '#FAF4D7',
    '#DED226',
    '#FCB370',
    '#ED8627',
    '#B65950',
    '#D96228',
    '#949935',
    '#628A29',
    '#1FA364',
    '#15ACA3',
    '#1E737C',
    '#0C4654',
    '#1E3542',
    '#0B2737',
    '#062012',
    '#1A0F28',
    '#0C0A0E',
]

for fn in argv[1:]:

    if not match(r'.*\.tif$', fn): continue

    gd = gdal.Open(fn)
    gt = gd.GetGeoTransform()
    nx = gd.RasterXSize
    ny = gd.RasterYSize
    arr = gd.ReadAsArray()
    hobs = arr[0,:,:]
    sky = arr[1:,:,:]

    fs = ( 8 * np.sqrt(nx/float(ny)), 6 / np.sqrt(nx/float(ny)) )

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
    m = ax.contourf(x, y, i2m(sky[3,:,:]), levels = clevels, colors = ccolors)
    plt.colorbar(m, ax = ax)
    ax.set_aspect('auto')

    plt.savefig(sub(r'\.tif$','.png',fn), dpi = 144)

    #---------------------------------------------------------------------------

    fig,ax = plt.subplots(figsize = fs)
    m = ax.contourf(x, y, i2m(sky[2,:,:]), levels = clevels, colors = ccolors)
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

    #---------------------------------------------------------------------------
