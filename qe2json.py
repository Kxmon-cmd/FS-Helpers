# This script takes in a single quantum espresso OUT-File and converts it into a json file(s) that can be read into fitsnap.
# The name of the Out-File and of the json file is set in the command line -> python qe2json.py -i inputfile.out -o output.json


import re
import json
import argparse
import numpy as np


def qe_out(filename, outname):
    """
    This funktion takes in the input file and the name of the desired output file. It will run throu the OUT-file and extrakt all the necissary information.
    The input of this funktion is a quantum espresso output file containing Energy, Force, Positions, ... . The output of this funktion is a dictionary containing
    all relevant information, which will be dumped into a json file.
    The positions and forces of all iterations are stored in an array and extracted later on.
    """

    NumAtoms = []
    Energy = []
    AtomTypes = []
    Positions = []
    Forces = []
    Lattice = []

    Forces_ = []
    Positions_ = []
    m=0

    with open(filename, "r") as f:
        lines = f.readlines()

    #start reading through every line in the input and search for the information

    for i, line in enumerate(lines):
        
        #reading the number of iterations done in the MD run
        if "nstep" in line:
            iterations = int(line[40:50])

        if "lattice parameter (alat)" in line:
            latpar = float(line[39:47].strip())*0.52917721090  #converting from bohr to angstom                        

        #reading amount of Atoms in the cell 
        if "     number of atoms/cell" in line:                     
            NumAtoms = [int(line[35:50].strip())]

        #reading atom types
        if "tau" in line:
            AtomTypes.append(str(line[17:21].strip()))

        #getting Lattice parameter
        if "               a" in line:
            a = [float(line[25:36].strip())*latpar, float(line[36:47].strip())*latpar, float(line[47:56].strip())*latpar]
            Lattice.append(a)

        #reading energies for every configuration and saving into an array 
        if "!    total energy" in line:                             
            Energy.append(float(line[35:50].strip())*13.6057039763)     #also converting rydberg to eV
        
        #reading forces for every atom
        if "1   force" in line:                                     
            f = [float(line[35:50].strip())*25.7110542481, float(line[50:63].strip())*25.7110542481, float(line[63:75].strip())*25.7110542481] #Converting Ry/au to eV/Angstrom
            Forces_.append(f)

        #reading positions of atoms
        #first extracting the different atom types
        Types = list(dict.fromkeys(AtomTypes))
        
        #searching for each type
        for t, type in enumerate(Types):

            if str(Types[t]) + "            "  in line:   
                #skipping first two lines                          
                m = m+1
                if m > 1:
                    Pos = [float(line[15:30].strip())*latpar, float(line[30:50].strip())*latpar, float(line[50:70].strip())*latpar]    #Converting from au to angstrom
                    Positions_.append(Pos)
    
    np.savetxt("Energy.txt",Energy[i],delimiter=',')

    #Extracting the Forces and Positions for one iteration
    for n in range(0,iterations):
        for j in range(0+n*NumAtoms[0], NumAtoms[0]+n*NumAtoms[0]):
            Forces.append(Forces_[j])                               

        for j in range(0+n*NumAtoms[0], NumAtoms[0]+n*NumAtoms[0]):
            Positions.append(Positions_[j])                         

        #crating a dictionary containg all information from the n-th interation
        data = {
            "Dataset": {
                "LatticeStyle": "angstrom",
                "EnergyStyle": "electronvolt",
                "AtomtypeStyle": "chemicalsymbol",
                "Positionstyle": "angstrom",
                "ForcesStyle": "electronvoltperangstrom",
                "Data": [{
                    "NumAtoms": NumAtoms[0],
                    "Lattice": Lattice,
                    "Energy": Energy[n],
                    "AtomTypes": AtomTypes,
                    "Positions": Positions,
                    "Forces": Forces,
                    }]
                }
            }
        
        #creating a jeson file for the n-th iteration
        with open(str(n+4000)+outname, "w") as f:                     
            json.dump(data, f, indent=2)

        #resetting the forces and positions
        Forces = []
        Positions = []

def main():
    parser = argparse.ArgumentParser(description='Convert QE MD output to FitSNAP JSON.')
    parser.add_argument('-i', '--input', required=True, help='QE output file')
    parser.add_argument('-o', '--output', required=True, help='Output JSON file')
    args = parser.parse_args()

    data = qe_out(args.input, args.output)

if __name__ == '__main__':
    main()
