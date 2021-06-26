# Relevant imports
import sys
import numpy as np
from Lattice import CahnHilliard as Lattice
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
# Figure Set up
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)

################################################################################

""" This program is for running the Cahn-Hilliard simulation for  varying
    inputs. It is able to run a free-energy collection mode [CH.Run()] or
    a parameter test for dx, dt [CH.TestParameters()]."""

################################################################################

if __name__=="__main__":

    # Handle command line inputs
    if len(sys.argv) == 1:
        print(f'Usage: python {sys.argv[0]} <Grid Size> <Phi-0> <StepNumber / Parameter Interval> <Animate>')
        quit()

    # Handle inputs
    GridSize = int(sys.argv[1])
    Phi0 = float(sys.argv[2])

    # Check Animation
    if sys.argv[3].upper() not in ['Y', 'N']:
        print('Invalid Input. Put in Y/N for Animation.')
    Animate = True if sys.argv[3].upper() == 'Y' else False

    # Generate Lattice
    CH = Lattice(GridSize, Phi0, 0.1, 0.1, 0.1, 0.1, Animate=True)

    # Which Test to perform
    TestType = int(input('1: Single Run\n2: Parameter Test\nSelect a test to run: '))
    # Free Energy
    if TestType == 1:

        # Run Simulation with Step Number
        FreeEnergyDensity, Steps = CH.Run(int(sys.argv[3]))

        # Plot Free Energy Figure
        plt.figure(figsize=(8,6))
        plt.xlabel('Steps')
        plt.ylabel('Specific Free Energy')
        plt.title('Plot of Specific Free Energy ')
        plt.plot(Steps, FreeEnergyDensity, color='black', lw=.5)
        plt.savefig(f'Plots/Cahn-Hilliard/{CH.Label}/SpecificFreeEnergyDensity[Phi0_{CH.Phi0}]_t{StepNumber}.png')
        plt.show()

        # Record Data
        OutputTable = Table([Steps, FreeEnergyDensity], names=['step', 'f'])
        OutputTable.write(f'./Data/Cahn-Hilliard/FreeEnergyDensity[Phi0_{Phi0}]_t{StepNumber}.dat', format='ascii', overwrite=True)

    # Parameter Test
    elif TestType == 2:

        # Run Simulation with Interval
        dx, dt, steps = CH.TestParameters(float(sys.argv[3]))
        # Record Data
        OutputTable = Table([dx, dt, steps], names=['dx', 'dt', 's'])
        OutputTable.write(f'./Data/Cahn-Hilliard/ParameterData/Phi0[{CH.Phi0}].dat', format='ascii', overwrite=True)

    # Error handling
    else:
        print('Invalid Input.')
        quit()
