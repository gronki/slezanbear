
import matplotlib.pyplot as plt
from numpy import pi, linspace, meshgrid, sin, cos, sqrt
from threshold import th as thr

th = linspace(0, 2 * pi, 1024)

x0 = linspace(-2, 2, 2**8)
y0 = linspace(-2, 2, 2**8)
x,y = meshgrid(x0, y0)

g = lambda costh: 0.75 * (1 + costh**2)
from numpy import where

muh, mul = 0.5, -0.5

print mul**2 * (2 * mul - 3 * muh) / (muh - mul) / 6 + (mul**2 - 1) / 2

lamp_jnorm = (2 + mul + muh) / 4
lamp_flxdn = (3 * (muh - mul) + mul**3) / (6 * (muh - mul) * lamp_jnorm)
lamp_flxup = muh**3 / (6 * (muh - mul) * lamp_jnorm)
print lamp_flxdn, lamp_flxup
glamp = lambda costh: where(costh < mul, 1, where(costh > muh, 0, (muh - costh) / (muh - mul))) / lamp_jnorm

refl_flxup = 4.0 / 3.0
grefl = lambda costh: 4 * where(costh > 0, costh, 0)
guni = lambda costh: 1.0 + 0*costh

integrate = lambda y,x: sum( (y[1:] + y[:-1]) * (x[1:] - x[:-1]) ) / 2

xx = cos(linspace(pi,0,2**12))
xxu = cos(linspace(pi/2,0,2**12))
xxd = cos(linspace(pi,pi/2,2**12))

print "Fu(lamp) = {:.5f}".format(integrate(glamp(xxu)*xxu,xxu))
print "Fu(refl) = {:.5f}".format(integrate(grefl(xxu)*xxu,xxu))

def fluxu(g):
    return "J = {:.5f}, F = {:.5f}, Fu = {:.5f}, Fd = {:.5f}".format(
            0.5 * integrate(g(xx),xx),
            integrate(g(xx)*xx,xx),
            integrate(g(xxu)*xxu,xxu),
            integrate(g(xxd)*xxd,xxd)
        )

g0 = lambda costh: 1.0 * where(costh > 0, 1, 0)
g1 = lambda costh: 1.5 * costh * where(costh > 0, 1, 0)
g2 = lambda costh: 6 * costh * (1 - costh) * where(costh > 0, 1, 0)
g3 = lambda costh: 15 * costh * (1 - costh)**2 * where(costh > 0, 1, 0)
g4 = lambda costh: 30 * costh * (1 - costh)**3 * where(costh > 0, 1, 0)

alb = 0.15
g = lambda costh: 0.5 * (glamp(costh) + alb * (lamp_flxdn / refl_flxup) * grefl(costh)) / (lamp_flxup + alb * lamp_flxdn)

a = 0.3
g = lambda costh: (1 - a) * g1(costh) + a * g4(costh)

print "LAMP " + fluxu(glamp)
print "REFL " + fluxu(grefl)
print "UNIF " + fluxu(guni)
print "L+R " + fluxu(g)


fig, ax = plt.subplots()
p = ax.pcolor(x0, y0, g(y / sqrt(x**2 + y**2)), cmap = 'hot', vmin = 0)
plt.colorbar(p)
ax.plot(g(cos(th)) * sin(th), g(cos(th)) * cos(th), color = '#053326')
ax.set_aspect(1)
plt.show()
