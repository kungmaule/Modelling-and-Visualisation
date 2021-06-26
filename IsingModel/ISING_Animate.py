# Import relevant libraries
import numpy as np
import sys, random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from astropy.table import Table

# Imports from Classes
from Lattice import Lattice

##### DESCRIPTION ##############################################################

""" This version of ISING.py is simply for animation purposes, methods for
    calculating values have been left in, but are not being stored, simply
    output to the terminal. Any other methods will be removed.             """

##### MAIN METHOD ##############################################################

if __name__=="__main__":

    # Check inputs for correct parameters
    if len(sys.argv) != 5:
        print(f'Incorrect Usage. \nCorrect Usage: python {sys.argv[0]} <G/K> <N> <Steps> <T>')
        quit()

    # Take in user inputs
    if sys.argv[1].upper() == 'G':
        Method = 'Glauber'
    elif sys.argv[1].upper() == 'K':
        Method = 'Kawasaki'
    else:
        print('Invalid Method Input. \nPlase input "K" for Kawasaki or "G" for Glauber Dynamics.')
        quit()
    GridLength = int(sys.argv[2])
    StepNumber = int(sys.argv[3])
    Temperature = float(sys.argv[4])

    # Hard-coded Parameters
    J = 1.0
    N = GridLength * GridLength
    EqSteps = 100
    Tau = 10

    # Define the Lattice
    print('Initialising Lattice...')
    IsingLattice = Lattice(GridLength, Method, J, Temperature)

    # Initialise animation
    fig = plt.figure(figsize=(6,6))
    col = cm.get_cmap('viridis', 4)
    newcol = ListedColormap(col(np.linspace(0.4,0.6,256)))
    im = plt.imshow(IsingLattice.Spins, animated=True, cmap=newcol)

    # Assign Temperature and update lattice
    print(f'Current Temperature: {Temperature}')

    # Iterate over Number of Steps
    for n in range(StepNumber):

        # Perform Sweep
        IsingLattice.Sweep()

        # Allow 100 steps to equilibriate
        if n <= EqSteps:

            # Update Plot every Tau Steps
            if (n%Tau) == 0:
                plt.cla()
                im=plt.imshow(IsingLattice.Spins, animated=True,cmap=newcol)
                plt.draw()
                plt.pause(0.0001)
                print(f'Step: {n}   Equilibrating')

            # Then continue
            continue

        # Otherwise record data every 10th step
        if (n%Tau) == 0:

            plt.cla()
            im=plt.imshow(IsingLattice.Spins, animated=True,cmap=newcol)
            plt.draw()
            plt.pause(0.0001)

            # Determine Lattice Energy
            Energy = IsingLattice.DetermineEnergy()

            # Check Method
            if IsingLattice.Method != 'Kawasaki':

                # Only Determine Magnetisation in Glauber Regime
                Magnetisation = IsingLattice.DetermineMagnetisation()

                print(f'Step: {n}   Energy: {Energy}    Magnetisation: {Magnetisation}')

            else:
                print(f'Step: {n}   Energy: {Energy}')

    # Finish Simulation
    print('Simulation Completed.')
