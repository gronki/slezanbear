# coding: utf-8

i1 = 2
i2 = 1
iqual = 6

from numpy import argsort, loadtxt, percentile, sqrt, exp, array

d = loadtxt('calib.txt', skiprows = 1)
s = argsort(d[:,iqual])

labels = ['alpha', 'relabs', 'hscale', 'beta', 'skybg']

import matplotlib.pyplot as plt

s1 = filter(lambda i: d[:,iqual][i] < percentile(d[:,iqual], 10), s)
s2 = filter(lambda i: d[:,iqual][i] < percentile(d[:,iqual], 2), s)
s3 = filter(lambda i: d[:,iqual][i] < percentile(d[:,iqual], 1), s)

plt.figure(figsize = (13.5,9))

for i, i1, i2 in [
    (1, 0, 2),
    (2, 1, 2),
    (3, 0, 1),
    (4, 0, 4),
    (5, 1, 4),
    (6, 0, 3),
]:
    plt.subplot(2, 3, i)
    plt.plot(d[:, i1], d[:, i2], '.', color = '#E5EBED')
    plt.plot(d[s1,i1], d[s1,i2], '.', color = '#E3C214')
    plt.plot(d[s2,i1], d[s2,i2], '.', color = '#EE4721')
    plt.plot(d[s3,i1], d[s3,i2], '.', color = '#300900')
    plt.xlabel(labels[i1])
    plt.ylabel(labels[i2])

plt.show()
