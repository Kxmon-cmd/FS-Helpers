import numpy as np
from numpy.polynomial.chebyshev import Chebyshev

#todo: add support for SBessel

#calculate cheb Value for argument x
def Cheb(n, x):
    if n == 0:
        Tn = 1.0
        cheb_val = Tn
    else :
        Tn = Chebyshev.basis(n)
        cheb_val = 0.5 - (0.5 * Tn(x))
    return cheb_val

def ChebExpCos(coeffs, r, rcut, rcut_in, dcut_in, lam):
    gr = []

    #set argument
    x =  1 - (2 * ((np.exp(-lam * (r / rcut)) - np.exp(-lam)) / (1 - np.exp(-lam))))
    for z, c in enumerate(coeffs):
        #c[1] = coeff, c[0] = n
        gr.append(c[1]*Cheb(c[0], x))
    return gr

def ChebPow(coeffs, r, lam, rcut):
    gr = []

    x = 2 * (1 - (1 - r/rcut)**lam) - 1
    for z, c in enumerate(coeffs):
        gr.append(c[1]*Cheb(c[0], x))
    return gr

def ChebLinear(coeffs, r, lam, rcut):
    gr = []

    x = (1 - r/rcut)
    for z, c in enumerate(coeffs):
        gr.append(c[1]*Cheb(c[0], x))
    return gr

def sinc(x):
    #define sinc function
    return np.sin(x) / x

def fn(n, r, rcut):
    #auxilliary function for SBessel
    return ((-1)**n * np.sqrt(2) * np.pi / rcut**1.5 *
            (n + 1) * (n + 2) / np.sqrt((n + 1)**2 + (n + 2)**2) *
            (sinc(r * (n + 1) / rcut) + sinc(r * (n + 2) / rcut)))


def SBessel(coeffs, rcut, r):
    """
    Simplified Bessel Basis (nur gr, gewichtet mit coeffs).
    coeffs: Liste von [n, faktor]
    """
    gr = []

    if r < rcut and r > 0:
        gr.append(fn(0, r, rcut))

        d_prev = 1.0
        for i in range(1, len(coeffs)):
            n, f = coeffs[i]
            en = n**2 * (n + 2)**2 / (4 * (n + 1)**4 + 1)
            dn = 1 - en / d_prev

            val = 1/np.sqrt(dn) * (fn(n, r, rcut) + np.sqrt(en/d_prev) * gr[i-1])  

            gr.append(val)
            d_prev = dn

    for m in range(len(gr)):
        gr[m] = gr[m] * coeffs[m][1]

    return gr