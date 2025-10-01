import numpy as np
from bin.cutoff import inner, cutoff, env
from bin.basis import ChebExpCos, ChebPow, ChebLinear, SBessel

def bit2int(a):
    #calculates the integer of a binary list
    return int(''.join(map(str, a)), 2) #+ 1

def ModCoeff(coeffs, basemax, a):
    #modify coeffs to drop irrelevant information

    coeff = [[idx, values[0]] for idx, values in coeffs]

    if bit2int(a) == 0:
        return coeff[0:basemax[bit2int(a)]]
    else:
        return coeff[np.sum(basemax[:bit2int(a)]):np.sum(basemax[:bit2int(a)]) + basemax[bit2int(a)]]
    

def V_r(r, rcut, dcut, coeffs, rcut_in, dcut_in, lam, base):
    #calculates the potential function
    fc = 1 - inner(rcut_in, dcut_in, r)

    #use corresponding basis function
    if base[0] == True:
        gr = ChebExpCos(coeffs, r, rcut, rcut_in, dcut_in, lam)
    if base[1] == True:
        gr = ChebPow(coeffs, r, rcut, rcut_in, dcut_in, lam)
    if base[2] == True:
        gr = ChebLinear(coeffs, r, lam, rcut)
    if base[3] == True:
        gr = SBessel(coeffs, rcut, r)

    if base[0] == True:
        #calculates the potential function with the coeffs and the basis functions. A Cutoff is also applied
        if r > (rcut - dcut) and r <= rcut:
            return np.sum(gr) * env(r, rcut) * cutoff(r, rcut, dcut)
        if r > rcut:
            return 0
        else:
            return np.sum(gr) * env(r, rcut) * fc
    else:
        return np.sum(gr) * fc