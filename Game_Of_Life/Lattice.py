# Import relevant libraries
import numpy as np
import sys, random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from astropy.table import Table

##### DESCRIPTION ##############################################################

""" This file contains the code for the Lattice Object, which is used to run the
    Game of Life simulations. """

##### CLASSES ####################################################################

class Lattice(object):

    """ This class contains all the attributes and methods for a Game of Life board,
    including standard running methods and set up methods. """

    ##### SET-UP METHODS #######################################################

    def __init__(self, GridLength, State, Animate=True):

        """ This method initialises the lattice with given parameters, and
        generates an intiail random spin layout. """

        # Set Parameters
        self.GridLength = GridLength
        self.Animate = Animate

        # Intialise Lattice
        self.Grid = np.zeros(shape=(GridLength, GridLength), dtype=float)
        # If Initial is None, initialise states randomly.
        if State == 'Randomised':
            self.RandomiseLattice()
        else:
            self.InitialiseState(State)

        self.SetInitialVariables()
        if self.Animate == True:
            self.InitaliseAnimation()

    def RandomiseLattice(self):
        """ This method randomly populates the lattice with alive and dead cells. """
        self.Grid = np.random.choice([0,1], size=(self.GridLength, self.GridLength))

    def InitialiseState(self, State):
        """ This method initialises the lattice in a specific formation according to
        a hardcoded state in a .dat file. """
        Start = int(self.GridLength/2.0)
        Init = Table.read(f'PreCodedStates/{State}.dat', format='ascii')
        iCoord = Init['x']
        jCoord = Init['y']
        for loc in range(iCoord.shape[0]):
            # Lobster and ACP144 are in a different format, working off absolutes, so are read in different.
            if State == 'Lobster' or State == 'ACP144':
                self.Grid[iCoord[loc]][jCoord[loc]] = 1.0
            else:
                self.Grid[Start+iCoord[loc]][Start+jCoord[loc]] = 1.0

    def SetInitialVariables(self):
        """ This method clears all current lattice variables and sets them to initial. """
        self.Count = np.zeros(shape=(self.GridLength, self.GridLength), dtype=float)
        self.GridNew = np.zeros(shape=(self.GridLength, self.GridLength), dtype=float)
        self.Stable = 0
        self.Oscillating = 0

    def Reset(self):
        """ This method resets all the variables of the lattice and resets the lattice. """
        self.SetInitialVariables()
        self.RandomiseLattice()

    def InitaliseAnimation(self):
        """ This method initialises the animation and colourbar. """
        fig = plt.figure(figsize=(6,6))
        plt.axis('off')
        col = cm.get_cmap('spring', 4)
        self.cmap = ListedColormap(col(np.linspace(0.3,0.9,256)))
        im = plt.imshow(self.Grid, animated=True, cmap=self.cmap, origin='lower', vmin=0, vmax=1)
        plt.colorbar(ticks=[0,1],fraction=0.046, pad=0.04)

    def AnimateFrame(self):
        """ This method updates the output frame. """
        plt.cla()
        plt.axis('off')
        im=plt.imshow(self.Grid, animated=True,cmap=self.cmap, origin='lower')
        plt.draw()
        plt.pause(0.0001)

    ##### METHODS ##############################################################

    def CountNeighbours(self, i, j):
        """ This method counts the number of nearest neighbours to a cell, including
        diagonals, and returns the number of them. """
        # Count Number of Alive Variables
        Count = 0
        for x in range(i-1, i+2):
            for y in range(j-1, j+2):
                # The count excludes the cell counting itself
                if (x,y) != (i,j):
                    Count += self.Grid[self.pbc(x, 0), self.pbc(y, 1)]
        return Count

    def Fate(self, i, j):
        """ This method determines the fate of a single cell. It counts the number
        of adjacent cells which are alive, and updates the new grid accordingly."""
        CurrentState = self.Grid[i,j]
        Count = self.CountNeighbours(i,j)

        if CurrentState == 0: # Dead
            if Count == 3:
                self.GridNew[i,j] = 1
            else:
                self.GridNew[i,j] = 0
        else: # Live
            if 2 <= Count <= 3:
                self.GridNew[i,j] = 1
            else:
                self.GridNew[i,j] = 0

    def Sweep(self):
        """ This method performs a sweep on the lattice. It determines the Fate
        of each individual lattice cell, and updates the grid. It also counts the
        number of alive cells before and after the sweep, and compares. From this
        it updates the stability or oscillatory variables accordingly.  """

        # Count number of alive cells.
        Alive = self.CountAlive()

        # Iterate over all cells
        for i in range(self.GridNew.shape[0]):
            for j in range(self.GridNew.shape[1]):
                self.Fate(i,j)

        # Reset the Lattice Grid (with new grid)
        self.ResetGrid()

        # Count number of living cells
        NewAlive = self.CountAlive()

        # Check for Oscillators of up to 2-site changes
        if 0 < abs(Alive-NewAlive) <= 2:
            self.Oscillating += 1
            self.Stable = 0
        # Check for Stability
        elif abs(Alive-NewAlive) == 0:
            self.Oscillating = 0
            self.Stable += 1
        # Reset State Variables
        else:
            self.Stable = 0
            self.Oscillating = 0

        # Animate if True
        if self.Animate == True:
            self.AnimateFrame()

    def pbc(self, x, axis):
        """ This method utilises periodic boundary conditions for a specific axis. """
        return np.mod(x, self.Grid.shape[axis])

    def CountAlive(self):
        return np.sum(self.Grid)

    def ResetGrid(self):
        """ This method sets the old grid to the new grid, and resets the new grid. """
        self.Grid = np.copy(self.GridNew)
        self.GridNew = np.zeros(shape=(self.GridLength, self.GridLength))

    def CenterOfMass(self):
        """ This method determines the center of mass of any alive cells on the
        grid, e.g. a spaceship, and returns it. """
        # Define Arrays
        x, y = [], []

        # Iterate over all cells
        for i in range(self.GridLength):
            for j in range(self.GridLength):

                # Check if Cell is alive
                if self.Grid[i][j] == 1.0:

                    # Check to see if cell is on the boundary
                    if (i == 0) or (i == self.GridLength-1):
                        return False
                    elif (j == 0) or (j == self.GridLength-1):
                        return False

                    # Otherwise record positon
                    else:
                        x.append(i)
                        y.append(j)

        return np.array([np.mean(x), np.mean(y)])

    ##### RUNNING METHODS ######################################################

    def RandomRun(self, DataPoints):

        """ This method runs a random lattice until it reaches stability a set number
        of times (defined by DataPoints) and records the data to Data/HistogramData.txt """

        # Hardcoded Values
        OscillationLimit = 100
        StabilityLimit = 10
        MaxSweeps = 5000
        count = 0

        # Iterate until the required number of data points have been collected
        while count < DataPoints:

            # Iterate over the maximum number of sweeps
            for sweep in range(MaxSweeps):

                # Sweep Lattice
                self.Sweep()

                # Output current lattice conditions every ten sweeps
                if sweep % 10 == 0:
                    print(f'{sweep}     Stable: {self.Stable}      Oscillating: {self.Oscillating}', end='\r')

                # Check for Stability
                if self.Stable == StabilityLimit:
                    print('\nStable at: ', sweep-10, '\n')
                    count += 1
                    # Write out Point
                    with open('Data/HistogramData.txt', 'a') as file:
                        file.write(f'{sweep-10}\n')
                    # Reset lattice + Break out of Loop
                    self.Reset()
                    break

                # Check for Oscillating behaviour
                elif self.Oscillating == OscillationLimit:
                    print(f'\nStuck in Oscillating State at: {sweep-200}')
                    # Reset Lattice + break out of loop
                    self.Reset()
                    break

                # Check for final state
                elif sweep == MaxSweeps-1:
                    print(f'Lattice did not stabilise within {MaxSweeps} sweeps.')
                    self.Reset()

    def OscillatorRun(self, NumberOfSteps):
        """ This method simply visualises an oscillator state. """
        # Iterate over number of steps
        for step in range(NumberOfSteps):
            self.Sweep()

    def SpaceshipRun(self, NumberOfSteps):
        """ This method visualises a spaceship state, and records all CoM data to
        a center of mass .dat file. """
        # Define Arrays
        CentreOfMass = np.zeros(shape=(0,2))
        Time = np.zeros(0)

        # Iterate over Number of Steps
        for step in range(NumberOfSteps):

            # Sweep Lattice
            self.Sweep()

            # Check for boundary crossing
            if np.any(self.CenterOfMass() == False):
                print('                                                             ', end='\r')
                print(f'Center of Mass: Crossing Boundary...', end='\r')
            else:
                CurrentCoM = self.CenterOfMass()
                CentreOfMass = np.vstack((CentreOfMass, CurrentCoM))
                Time = np.append(Time, step)
                print('                                                             ', end='\r')
                print(f'Center of Mass: ({CurrentCoM[0]},{CurrentCoM[1]})', end='\r')

            # Write out Data
        CoMData = Table([Time, CoM[:,0], CoM[:,1]], names=['Time', 'X', 'Y'])
        CoMData.write(f'Data/CenterOfMass({State}).dat', format='ascii', overwrite=True)
