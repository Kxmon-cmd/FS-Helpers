import re
import numpy as np
import matplotlib.pyplot as plt
import argparse
from bin.basis import ChebExpCos, ChebPow, ChebLinear, SBessel
from bin.Pot import V_r, ModCoeff, bit2int

def main():

    #Adding Pot file as input
    parser = argparse.ArgumentParser(description='Plot two-body term')
    parser.add_argument('-i', '--input', required=True, help='ACE Potential file')
    args = parser.parse_args()

    Pot = args.input

    with open(Pot, "r") as f:
        lines = f.readlines()

    #bond is defined here
    a = [0, 0] 

    base = []
    coeffs = []
    rcut = []
    dcut = []
    dcut_in = []
    rcut_in = []
    lam = []
    basemax = []

    #extracting information
    for i, line in enumerate(lines):            
        if "rank: 1" in line and "ls: [0]" in line:
                    ct_match = re.search(r'ctildes: \s*\[(.*?)\]', line)
                    ns_match = re.search(r'ns: \[(\d+)\]', line)
                    ct = ct_match.group(1)
                    ct = [float(x.strip()) for x in ct.split(',')]
                    n = int(ns_match.group(1))
                    coeffs.append([n, ct])

        if "elements: " in line:
            match = re.search(r'elements:\s*\[(.*?)\]', line)
            ele = match.group(1)
            elements = [x.strip() for x in ele.split(',')]

        if "nradbasemax: " in line:
            match = re.search(r'nradbasemax:\s*([0-9.eE+-]+)', line)
            basemax.append(int(match.group(1)))

        if "rcut: " in line:
            rcut_match = re.search(r'\s*rcut:\s*([0-9.eE+-]+)', line)
            rcut.append(float(rcut_match.group(1)))

        if "dcut: " in line:
            dcut_match = re.search(r'\s*dcut:\s*([0-9.eE+-]+)', line)
            dcut.append(float(dcut_match.group(1)))

        if "rcut_in: " in line:
            rcut_in_match = re.search(r'\s*rcut_in:\s*([0-9.eE+-]+)', line)
            rcut_in.append(float(rcut_in_match.group(1)))

        if "dcut_in: " in line:
            dcut_in_match = re.search(r'\s*dcut_in:\s*([0-9.eE+-]+)', line)
            dcut_in.append(float(dcut_in_match.group(1)))
            
        if "lambdahc: " in line:
            lam_match = re.search(r'\s*lambdahc:\s*([0-9.eE+-]+)', line)
            lam.append(float(lam_match.group(1)))

        if "radbasename: " and str(a) + ":" in line:
            #check which basis is used
            CEC_match = re.search(r'\s*radbasename: ChebExpCos', line)
            CP_match = re.search(r'\s*radbasename: ChebPow', line)
            CL_match = re.search(r'\s*radbasename: ChebLinear', line)
            SB_match = re.search(r'\s*radbasename: SBessel', line)
            
            base = [False, False, False, False] #[CEC, CP, CL, SB]
            if CEC_match:
                base[0] = True
            if CP_match:
                base[1] = True
            if CL_match:
                base[2] = True
            if SB_match:
                base[3] = True
    
    coeff = ModCoeff(coeffs, basemax, a)        #modified coeff array according to bond defined in a
    i = bit2int(a)                              #number corresponding to bond defined in a 

    #plotting potential
    r_vals = np.linspace(0 , rcut[i], 3000)
    V_vals = [V_r(r, rcut[i], dcut[i], coeff, rcut_in[i], dcut_in[i], lam[i], base)/4 for r in r_vals]

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
    
