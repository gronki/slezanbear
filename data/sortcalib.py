# coding: utf-8

from numpy import argsort, loadtxt

d = loadtxt('calib.txt', skiprows = 1)
s = argsort(d[:,5])

for i in range(d.shape[0]):
    print "{d[0]:7.3f} {d[1]:6.3f} {d[2]:7.1f} {d[3]:5.3f} {d[4]:6.3f} {d[5]:8.4f} {d[6]:8.4f} {d[7]:8.4f}".format(d = d[s[i],:])
