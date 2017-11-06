from sympy import *

def getnorm(f):
    from sympy import integrate, symbols
    x = symbols('x')
    I = integrate(2 * f(x) * x, x)
    norm = 1 / (I.subs(x,1) - I.subs(x,0))
    return norm

print "{}".format(getnorm(lambda x: 1))

print "{} * x".format(getnorm(lambda x: x))
print "{} * x * (1 - x)".format(getnorm(lambda x: x * (1 - x)))
print "{} * x * (1 - x)**2".format(getnorm(lambda x: x * (1 - x)**2))
print "{} * x * (1 - x)**3".format(getnorm(lambda x: x * (1 - x)**3))
print "{} * x * (1 - x)**4".format(getnorm(lambda x: x * (1 - x)**4))
print "{} * x * (1 - x)**5".format(getnorm(lambda x: x * (1 - x)**5))
print "{} * x * (1 - x)**6".format(getnorm(lambda x: x * (1 - x)**6))
print "{} * x * (1 - x)**7".format(getnorm(lambda x: x * (1 - x)**7))

print "{} * (1 - x)".format(getnorm(lambda x: 1 - x))
print "{} * (1 - x)**2".format(getnorm(lambda x: (1 - x)**2))
print "{} * (1 - x)**3".format(getnorm(lambda x: (1 - x)**3))
print "{} * (1 - x)**4".format(getnorm(lambda x: (1 - x)**4))
print "{} * (1 - x)**5".format(getnorm(lambda x: (1 - x)**5))
print "{} * (1 - x)**6".format(getnorm(lambda x: (1 - x)**6))

print "{} * sqrt(x)".format(getnorm(lambda x: sqrt(x)))
print "{} * sqrt(x) * (1 - x)".format(getnorm(lambda x: sqrt(x) * (1 - x)))
print "{} * sqrt(x) * (1 - x)**2".format(getnorm(lambda x: sqrt(x) * (1 - x)**2))
print "{} * sqrt(x) * (1 - x)**3".format(getnorm(lambda x: sqrt(x) * (1 - x)**3))
print "{} * sqrt(x) * (1 - x)**4".format(getnorm(lambda x: sqrt(x) * (1 - x)**4))
