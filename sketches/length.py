from numpy import linspace, logspace
from numpy import sin, cos, log, sqrt, arccos, exp, abs
import matplotlib.pyplot as plt
import latexrc

#-------------------------------------------------------------------------------

def main(R,H,ax):

    Dmax = R * arccos(R / (R + H))
    D = linspace(0, Dmax, 2**10)

    L0 = sqrt(H**2 + D**2)
    L1 = sqrt(H**2 + (1 + H/R) * D**2)
    L2 = sqrt(H**2 + 2*R*(R+H)*(1 - cos(D/R)))

    #---------------------------------------------------------------------------

    ax.plot(D / R, L0 / H, color = '#B9BD00', label = 'flat earth')
    ax.plot(D / R, L1 / H, color = '#00A66F', label = 'round earth (approx)')
    ax.plot(D / R, L2 / H, color = '#131313', linestyle = '--', label = 'round earth (exact)')
    ax.set_title('$H / R = 1 / {}$'.format(int(R/H)))
    ax.set_xlim(0, Dmax / R)
    ax.set_ylim(1, None)
    ax.set_xlabel('$D / R$')
    ax.set_ylabel('$L / H$')

#-------------------------------------------------------------------------------

R = 6371.0
fig, axes = plt.subplots(1, 2, figsize = (6,3.2))
main(R, R / 2, axes[0])
main(R, R / 20, axes[1])
axes[1].legend()
plt.tight_layout()
plt.savefig('/tmp/length.pdf')
