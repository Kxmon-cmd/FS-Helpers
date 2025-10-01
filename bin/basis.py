import numpy as np
from numpy.polynomial.chebyshev import Chebyshev

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

def fn(n, r, rcut):
    return (-1)**n * np.sqrt(2) * np.pi / rcut**1.5  * (n + 2) / np.sqrt(np.sqrt(n + 1) + np.sqrt(n +2)) * (np.sinc(r * (n + 1) * np.pi / rcut) + np.sinc(r * (n + 2) * np.pi / rcut))

def SBessel(coeffs, rcut, r):
    gr = []

    if r > rcut:
        gr.append(fn(0, r, rcut))

        for z, c in enumerate(coeffs):
            if z > 0:
                en = np.sqrt(c[0]) * np.sqrt(c[0] + 2) / (4 * (c[0] + 1)**4 +1)
                dn = 1 - en

                x = 1/ np.sqrt(dn) * (fn(r, rcut, c[0]) + np.sqrt(en) * gr[c[0] - 1])
                gr.append()

    else:
        gr.append(0)
    
    return gr