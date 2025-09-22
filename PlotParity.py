import argparse
import matplotlib.pyplot as plt
from fitsnap3lib.tools.dataframe_tools import DataframeTools
import numpy as np

def main():

    #Adding df file as input
    parser = argparse.ArgumentParser(description='Plot Energy/Force Parity')
    parser.add_argument('-i', '--input', required=True, help="FitSNAP Dataframe")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-e", "--energy",
        action="store_true",
        help="Plot energy parity"
    )
    group.add_argument(
        "-f", "--force",
        action="store_true",
        help="Plot force parity"
    )

    args = parser.parse_args()
    df = args.input

    df = DataframeTools(df)
    df.read_dataframe()

    if args.energy:
        dfe  = df.df[df.df['Row_Type']=='Energy']
        x = np.linspace(dfe["preds"].min(), dfe["preds"].max())
    if args.force:
        dfe  = df.df[df.df['Row_Type']=='Force']
        x = np.linspace(dfe["preds"].min(), dfe["preds"].max())

    y = x
    dfe_Test  = dfe[dfe["Testing"]==True]
    dfe_Training = dfe[dfe["Testing"]==False]

    plt.scatter(dfe_Training['preds'],dfe_Training['truths'], c="gold", s=4, alpha=1)
    plt.scatter(dfe_Test['preds'],dfe_Test['truths'], c="cadetblue", s=4, alpha=1)
    plt.plot(x,y, c="teal")
    plt.xlabel(r"$E_{ML} (eV)$")
    plt.ylabel(r"$E_{ref} (eV)$")
    plt.show()

if __name__ == "__main__":
    main()