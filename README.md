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

