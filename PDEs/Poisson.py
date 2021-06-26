import sys
import numpy as np
from Lattice import Poisson as Lattice
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)

##### METHODS ##################################################################

def GetCoordinates(GridLength):

    """ Generates the monopole coordinates from user input. """

    Coords = input('Input Monopole Coordinates i j k: ')
    ijk = Coords.split(' ')
    i, j, k = ijk[0], ijk[1], ijk[2]
    if not 0 <= np.any([int(i),int(j),int(k)]) <= GridLength:
        print('Invalid Coordinates Given.')
        quit()
    return i,j,k

def AssumeCoordinates(GridLength):

    """ Assumes central monopole charge / wire"""

    return GridLength/2, GridLength/2, GridLength/2

##### MAIN #####################################################################

if __name__=="__main__":

    if len(sys.argv) == 1:
        print(f'Usage: python {sys.argv[0]} <Grid Size> <Phi-0> <Mode> <Field (E/B)> <Animate (Y/N)>')
        quit()
    # Handle inputs
    GridSize = int(sys.argv[1])
    Phi0 = float(sys.argv[2])
    Noise = 0.1

    # Check Run Mode
    if sys.argv[3].upper() not in ['J', 'GS', 'SOR']:
        print('Invalid Mode. Select J (Jacobi), GS (Gauss-Seidel) or SOR (overrelaxation)')
        quit()
    else:
        Mode = sys.argv[3].upper()
    # Check field
    if sys.argv[4].upper() not in ['E', 'B']:
        print('Invalid Field. Select E (Electric) or B (Magnetic)')
        quit()
    else:
        Field = sys.argv[4].upper()
    # Check Animation
    if sys.argv[4].upper() not in ['Y', 'N']:
        print('Invalid Input. Put in Y/N for Animation.')
    Animate = True if sys.argv[4].upper() == 'Y' else False

    # Jacobi Mechanism
    if Mode == 'J':

        print('\nRunning Jacobi Algorithm...\n')

        # Get Monopole Coordinates
        i, j, k = AssumeCoordinates(GridSize)
        # Generate Lattice
        Jacobi = Lattice(GridSize, Phi0, Noise=0.1, i=int(i), j=int(j), k=int(k), Field=Field, Animate=Animate)
        # perform Run
        Jacobi.JacobiRun()

    # Gauss Seidel Mode
    elif Mode == 'GS':

        print('\nRunning Gauss-Seidel Algorithm...\n')

        # Get Monopole Coordinates
        i, j, k = AssumeCoordinates(GridSize)
        # Generate Lattice
        GaussSeidel = Lattice(GridSize, Phi0, Noise=0.1, i=int(i), j=int(j), k=int(k), Field=Field, Animate=Animate)
        # Perform Run
        GaussSeidel.GaussSeidelRun()

    # SOR Mode
    elif Mode == 'SOR':

        print('\nRunning SOR Algorithm...\n')

        # Get Monopole Coordinates
        i, j, k = AssumeCoordinates(GridSize)
        # Generate Lattice
        SOR = Lattice(GridSize, Phi0, Noise=0.1, i=int(i), j=int(j), k=int(k), Field=Field, Animate=Animate)

        # Set up different test conditions
        TestType = int(input('1: Single ω\n2: Convergence\nWhich Test: '))

        # Single Run Test
        if TestType == 1:

            # Take input Omega and run
            omega = float(input('Input ω: '))
            SOR.SOR_Single(omega)

        # Run Convergence Test
        elif TestType == 2:

            # Take input Omega and run
            Interval = float(input('Input ω Interval: '))
            SOR.SOR_Convergence(Interval)
        # Error Handling
        else:
            print('Invalid test input.')
            quit()
