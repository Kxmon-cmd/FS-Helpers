import re
import numpy as np
import matplotlib.pyplot as plt
import argparse
from bin.cutoff import inner, cutoff, env
from bin.basis import ChebExpCos, ChebPow, ChebLinear

#Todo: add other basis functions

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
    

    #calculates the potential function with the coeffs and the basis functions. A Cutoff is also applied
    if r > (rcut - dcut) and r <= rcut:
        return np.sum(gr) * env(r, rcut) * cutoff(r, rcut, dcut)
    if r > rcut:
        return 0
    else:
        return np.sum(gr) * env(r, rcut) * fc
    
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
                    coeffs.append([n, c])


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

        if "radbasename: " in line:
            #check which basis is used
            CEC_match = re.search(r'\s*radbasename: ChebExpCos', line)
            CP_match = re.search(r'\s*radbasename: ChebPow', line)
            CL_match = re.search(r'\s*radbasename: ChebLinear', line)
            
            base = [False, False, False] #[CEC, CP, CL]
            if CEC_match:
                base[0] = True
            if CP_match:
                base[1] = True
            if CL_match:
                base[2] = True
    
    #plotting potential
    r_vals = np.linspace(0 , rcut, 3000)
    V_vals = [V_r(r, rcut, dcut, coeffs, rcut_in, dcut_in, lam, base)/4 for r in r_vals]

    plt.figure(figsize=(8, 5))
    plt.plot(r_vals, V_vals, label="Two-body term V(r)")
    plt.xlabel("r [Å]")
    plt.ylabel("V(r) [eV]")
    plt.title("ACE two-body term (ChebExpCos)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
    
