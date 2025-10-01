# FS-Helpers

This repository contains helpful Python scripts for working with [Quantum ESPRESSO](https://www.quantum-espresso.org/) and [FitSNAP](https://fitsnap.github.io/).  
The tools focus on converting MD outputs into FitSNAP-compatible formats and visualising results.

## Features
- **qe2json.py**: Convert MD output files into json format.
- **qe2xyz.py**: Convert MD output into xyz format.
- **PlotParity.py**: Plot predicted energies/forces vs. reference data from a `FITSNAP.df` dataframe.  
- **2bterm.py**: Visualise the two-body term of an ACE potential.    

## Usage

- **qe2json.py**:
  
       python qe2json.py -i qe.out -o outname.json
- **qe2xyz.py**:

      python qe2xyz.py -i qe.out -o outname.xyz
- **PlotParity.py**:
  
      python PlotParity.py -i FitSNAP.df -e
      python PlotParity.py -i FitSNAP.df -f
- **2bterm.py**:

      python 2bterm.py -i ACE_pot.yace

## Parameters

To use the "2bterm.py" on systems with multiple bonds, the parameter "a" can be used to change the bond which is being plotted. Note that the code currently only supports systems of up to two species.
Some example potentials can be accessed in the examples folder. (AlN Yang, G., Liu, Y.-B., Yang, L., & Cao, B.-Y. (2024). Machine-learned atomic cluster expansion potentials for fast and quantum-accurate thermal simulations of wurtzite AlN. Journal of Applied Physics, 135, 085105. https://doi.org/10.1063/5.0188905)

