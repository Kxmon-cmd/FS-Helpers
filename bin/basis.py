import numpy as np
from numpy.polynomial.chebyshev import Chebyshev

#todo: add the other basis functions 

def ChebExpCos(n, r, rcut, rcut_in, dcut_in, lam):

    #Calculates the chebyshev basis functions
    x =  1 - (2 * ((np.exp(-lam * (r / rcut)) - np.exp(-lam)) / (1 - np.exp(-lam))))
    if n == 0:
        Tn = 1
        cheb_val = Tn
    else :
        Tn = Chebyshev.basis(n)
        cheb_val = 0.5 - (0.5 * Tn(x))
    return cheb_val