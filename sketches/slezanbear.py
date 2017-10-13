# -*- coding: utf-8 -*-

from factory import *

__lib = CDLL('libslezanbear.so')

sbmap = CFACTORY(__lib, 'sbmap_c', [
    (ndpointer(c_float,2), 'in', 'mapi'),
    (c_int, 'in', 'nxi'),
    (c_int, 'in', 'nyi'),
    (ndpointer(c_float,1,6), 'in', 'gi'),
    (ndpointer(c_float,2), 'in', 'maph'),
    (c_int, 'in', 'nxh'),
    (c_int, 'in', 'nyh'),
    (ndpointer(c_float,1,6), 'in', 'gh'),
    (ndpointer(c_double,2), 'inout', 'I1'),
    (ndpointer(c_double,2), 'inout', 'I2'),
    (c_int, 'in', 'nx'),
    (c_int, 'in', 'ny'),
    (ndpointer(c_float,1,6), 'in', 'g'),
])
