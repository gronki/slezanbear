# coding: utf-8

from numpy import argsort, loadtxt, percentile

d = loadtxt('calib.txt', skiprows = 1)
s = argsort(d[:,6]**2 + d[:,7]**2)

i1 = 0
i2 = 2

labels = ['alpha', 'relabs', 'hscale', 'beta', 'skybg']

import matplotlib.pyplot as plt

s1 = filter(lambda i: d[i,5] < percentile(d[:,5], 10), s)
s2 = filter(lambda i: d[i,5] < percentile(d[:,5], 2), s)

plt.figure(figsize = (12,6))

plt.subplot(1, 2, 1)
plt.plot(d[:, i1], d[:, i2], '.', color = '#C8DBE2')
plt.plot(d[s1,i1], d[s1,i2], '.', color = '#C6A600')
plt.plot(d[s2,i1], d[s2,i2], '.', color = '#C62400')
plt.xlabel(labels[i1])
plt.ylabel(labels[i2])

R = [
    (min(d[:,0]), max(d[:,0])),
    (min(d[:,1]), max(d[:,1])),
    (min(d[:,2]), max(d[:,2])),
    (min(d[:,3]), max(d[:,3])),
    (min(d[:,4]), max(d[:,4])),
]

from numpy import linspace, meshgrid, sqrt
x1 = linspace(R[i1][0], R[i1][1], 48)
x2 = linspace(R[i2][0], R[i2][1], 48)
X1, X2 = meshgrid(x1, x2)
Y = X1 * 0

for i in range(x1.shape[0]):
    for j in range(x2.shape[0]):
        dist2 =   (d[:,i1] - X1[i,j])**2 / (R[i1][1] - R[i1][0])**2 \
                + (d[:,i2] - X2[i,j])**2 / (R[i2][1] - R[i2][0])**2
        # w = 1 / (dist2 + 0.01)
        w = 1 / (d[:,5] / min(d[:,5]) - 0.9) / (dist2 + 0.01)
        Y[i,j] = sum(w * d[:,5]) / sum(w)

plt.subplot(1, 2, 2)
plt.contourf(x1, x2, Y, 16, cmap = 'inferno_r')
plt.xlabel(labels[i1])
plt.ylabel(labels[i2])
plt.colorbar()

plt.show()
