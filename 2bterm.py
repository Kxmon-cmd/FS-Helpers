import re
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import Chebyshev
import argparse

def inner(rcut_in, dcut_in, r):

    #inner cutoff
    if r <= rcut_in - dcut_in:
        return 1
    elif r >= rcut_in:
        return 0
    else:
        x = 1 - 2 * (1 + (r - rcut_in) / dcut_in)
        return 0.5 + 7.5/2 * (x/4 - (x**3)/6 + (x**5)/20)

def phi_n(n, r, rcut, rcut_in, dcut_in, lam):

    #Calculates the chebyshev basis functions
    x =  1 - (2 * ((np.exp(-lam * (r / rcut)) - np.exp(-lam)) / (1 - np.exp(-lam))))
    if n == 0:
        Tn = 1
        cheb_val = Tn
    else :
        Tn = Chebyshev.basis(n)
        cheb_val = 0.5 - (0.5 * Tn(x))
    return cheb_val

def cutoff(r, rcut, dcut):
    #outer cutoff
    return 0.5 * (1.0 + np.cos(np.pi * (r - (rcut - dcut)) / dcut))

def V_r(r, rcut, dcut, coeffs, rcut_in, dcut_in, lam):
    #calculates the potential function

    fc = 1 - inner(rcut_in, dcut_in, r)

    #envelope function
    env = 0.5 * (1 + np.cos(np.pi * (r / rcut)))

    #calculates the potential function with the coeffs and the basis functions. A Cutoff is also applied
    if r > (rcut - dcut) and r <= rcut:
        return sum(c * phi_n(n, r, rcut, rcut_in, dcut_in, lam) for n, c in coeffs) * env * cutoff(r, rcut, dcut)
    if r > rcut:
        return 0
    else:
        return sum(c * phi_n(n, r, rcut,  rcut_in, dcut_in, lam) for n, c in coeffs) * env * fc
    
def main():

    #Adding Pot file as input
    parser = argparse.ArgumentParser(description='Plot two-body term')
    parser.add_argument('-i', '--input', required=True, help='ACE Potential file')
    args = parser.parse_args()

    Pot = args.input

    with open(Pot, "r") as f:
        lines = f.readlines()

    coeffs = []

    #extracting information
    for i, line in enumerate(lines):
        if "rank: 1" in line and "ls: [0]" in line:
                ns_match = re.search(r'ns: \[(\d+)\]', line)
                ct_match = re.search(r'ctildes: \[([-0-9.eE+]+)\]', line)
                if ns_match and ct_match:
                    n = int(ns_match.group(1))
                    c = float(ct_match.group(1))
                    coeffs.append((n, c))


        if "rcut: " in line:
            rcut_match = re.search(r'\s*rcut:\s*([0-9.eE+-]+)', line)
            rcut = float(rcut_match.group(1))

        if "dcut: " in line:
            dcut_match = re.search(r'\s*dcut:\s*([0-9.eE+-]+)', line)
            dcut = float(dcut_match.group(1))

        if "rcut_in: " in line:
            rcut_in_match = re.search(r'\s*rcut_in:\s*([0-9.eE+-]+)', line)
            rcut_in = float(rcut_in_match.group(1))

        if "dcut_in: " in line:
            dcut_in_match = re.search(r'\s*dcut_in:\s*([0-9.eE+-]+)', line)
            dcut_in = float(dcut_in_match.group(1))
            
        if "lambdahc: " in line:
            lam_match = re.search(r'\s*lambdahc:\s*([0-9.eE+-]+)', line)
            lam = float(lam_match.group(1))
            
    #plotting potential
    r_vals = np.linspace(0 , rcut, 3000)
    V_vals = [V_r(r, rcut, dcut, coeffs, rcut_in, dcut_in, lam)/4 for r in r_vals]

    plt.figure(figsize=(8, 5))
    plt.plot(r_vals, V_vals, label="Zweikörperpotential V(r)")
    plt.xlabel("r [Å]")
    plt.ylabel("V(r) [eV]")
    plt.title("ACE 2-Körperpotential (ChebExpCos)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
    
