from numpy import exp

def th(x, w = 1.0):
    return 1 / ( 1 + exp(-4*x/w) )
