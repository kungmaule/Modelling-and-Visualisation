# Import Relevant Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys, time
from Lattice import Lattice
from astropy.table import Table

##### DESCRIPTION ##############################################################

""" This method runs the simulation for a variety of input conditions. It can
    run a random lattice, or a pre-coded state dependent on user input.   """

##### MAIN METHOD ##############################################################

if __name__=="__main__":

    # Input Handling
    if len(sys.argv) != 4:
        print(f'Incorrect Usage.\nCorrect Usage: python {sys.argv[0]} <Size> <Randomise(Y/N)> <Animate(Y/N)>')
        quit()

    else:
        # Record Grid Length
        GridLength = int(sys.argv[1])

        # Check Simulation Mode
        if sys.argv[2].upper() == 'Y':
            State = 'Randomised'
        elif sys.argv[2].upper() == 'N':
            PreCodedStates = ['Glider', 'Blinker', 'Lwss', 'Penta', 'Pulsar', 'KokGalaxy', 'Lobster', 'ACP144']
            print(f'Select the number from the available states are:')
            for count, Name in enumerate(PreCodedStates):
                print(f'{count}: {Name}')
            n = None
            while n not in np.arange(0, len(PreCodedStates)):
                n = int(input('InitialState: '))
            State = PreCodedStates[n]
        else:
            print('Invalid Input. Put in Y/N for an initial randomised state.')
            quit()

        # Check Animation
        if sys.argv[3].upper() not in ['Y', 'N']:
            print('Invalid Input. Put in Y/N for Animation.')
        Animate = True if sys.argv[3].upper() == 'Y' else False

    # Initalise Lattice
    GameBoard = Lattice(GridLength, State, Animate=Animate)

    # CASE 1: Random Lattice for Stability Data Collection
    if State == 'Randomised':

        # Determine how many data points to collect
        DataPoints = int(input('Input number of data points to collect: '))

        # Run Simulation
        print('\nStarting Randomised Simulation...\n')
        GameBoard.RandomRun(DataPoints)
        print('\nSimulation Completed. \nStabilisation Steps are available under [Data/HistogramData.txt]')

        # Program Ends

    else:

        # Get user input number of steps
        StepNumber = int(input('Input number of steps to run: '))

        # If Spaceship, record CoM data, otherwise just visualise.
        if State in ['Glider', 'Lwss', 'Lobster']:
            GameBoard.SpaceshipRun(StepNumber)
        else:
            GameBoard.OscillatorRun(StepNumber)

        # Program Ends
