from ase.io import read, write
import os
import argparse

def qe2XYZ(filename, outname):

    # Read all frames from QE output
    try:
        traj = read(filename, index=":")
    except Exception as e:
        print("Error reading QE output file:", e)
        exit(1)

    # Check if energies and forces are present
    for i, atoms in enumerate(traj):
        if 'energy' not in atoms.info:
            print(f"[Warning] Frame {i} missing energy.")

        if 'forces' not in atoms.arrays:
            print(f"[Warning] Frame {i} missing forces.")

    #write xyz file
    write(outname, traj, format='extxyz')

def main():
    parser = argparse.ArgumentParser(description='Convert QE MD output to XYZ.')
    parser.add_argument('-i', '--input', required=True, help='QE output file')
    parser.add_argument('-o', '--output', required=True, help='Output XYZ file')
    args = parser.parse_args()

    data = qe2XYZ(args.input, args.output)

if __name__ == '__main__':
    main()
