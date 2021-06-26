# Import Relevant Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys, time
from Lattice import Lattice
from astropy.table import Table

##### DESCRIPTION ##############################################################

""" This method runs the SIRS simulation in various user defined modes. Acts as
    a 'main menu' for the code.                                              """

##### MAIN METHOD ##############################################################

if __name__=="__main__":

    RunModes = ['P', 'S', 'C', 'I']

    # Check User Input
    if len(sys.argv) != 3:
        print(f'Correct Usage: {sys.argv[0]} <N> <Mode>')
        #<S-I> <I-R> <R-S>')
        quit()
    else:
        GridLength = int(sys.argv[1])
        Mode = sys.argv[2]
        EqbSteps = 100
        N = GridLength**2
        if Mode.upper() not in RunModes:
            print('Avaialble Modes: \n<S> Single Run \n<P> Generate Plane Data \n<C> Generate p2,p3 Cut Data \n<I> Immune Fraction')
            quit()

    if Mode == 'S':
        StepNumber = int(input('Input Number of Steps: '))
        p1 = float(input('Input Probability of S -> I: '))
        p2 = float(input('Input Probability of I -> R: '))
        p3 = float(input('Input Probability of R -> S: '))

        # Check Animation
        AnimateInput = False
        while AnimateInput not in [1,2]:
            AnimateInput = int(input('Animate? (1: Yes, 2: No): '))
        Animate = True if AnimateInput == 1 else False

        # Run Single Run
        SIRS = Lattice(GridLength, Animate=Animate)
        SIRS.SingleRun(StepNumber, p1,p2,p3)

    elif Mode == 'P':
        StepNumber = int(input('Input Number of Steps: '))
        pInterval = float(input('Input probability interval for p1: '))
        Animate = False # Assume Animation False

        # Run
        SIRS = Lattice(GridLength, Animate=Animate)
        SIRS.Generate_p1p3_Plane(StepNumber, pInterval, p2=0.5)

    elif Mode == 'C':
        StepNumber = int(input('Input Number of Steps: '))
        pInterval = float(input('Input probability interval for p1: '))
        Animate = False # Assume Animation False

        # Run
        SIRS = Lattice(GridLength, Animate=Animate)
        SIRS.FindWaves(StepNumber, pInterval, 0.2, 0.5, 0.5, 0.5)

    elif Mode == 'I':
        StepNumber = int(input('Input Number of Steps: '))
        fImmInterval = float(input('Input interval for immune fraction: '))
        p1 = float(input('Input Probability of S -> I: '))
        p2 = float(input('Input Probability of I -> R: '))
        p3 = float(input('Input Probability of R -> S: '))
        Animate = False

        # Run
        SIRS = Lattice(GridLength, Animate=Animate, Immunity=True)
        SIRS.ImmunityRun(StepNumber, fImmInterval, p1, p2, p3)
