import numpy as np

#todo: Add support for zbl

def inner(rcut_in, dcut_in, r):
    #inner cutoff
    if r <= rcut_in - dcut_in:
        return 1
    elif r >= rcut_in:
        return 0
    else:
        x = 1 - 2 * (1 + (r - rcut_in) / dcut_in)
        return 0.5 + 7.5/2 * (x/4 - (x**3)/6 + (x**5)/20)
    
def cutoff(r, rcut, dcut):
    #outer cutoff
    return 0.5 * (1.0 + np.cos(np.pi * (r - (rcut - dcut)) / dcut))

def env(r, rcut):
    #envolope function
    return 0.5 * (1 + np.cos(np.pi * (r / rcut)))