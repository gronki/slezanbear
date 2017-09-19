from numpy import exp, where, abs

def erf(x):
    p  =  0.327591100
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    x0 = abs(x)
    t = 1 / (1 + p*x0)
    erf0 = 1 - t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5)))) * exp(-x0**2)
    return where(x >= 0, erf0, -erf0)

if __name__ == '__main__':
    from math import erf as matherf
    from numpy import array,linspace
    import matplotlib.pyplot as plt
    x = linspace(-4,4,500)
    fig, axes = plt.subplots(2,1,figsize=(6,10))
    axes[0].plot(x, erf(x), label = 'appx')
    erfm = array([matherf(xi) for xi in x])
    axes[0].plot(x, erfm, label = 'math')
    axes[0].legend()
    axes[1].plot(x, erf(x) - erfm)
    plt.show()
