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