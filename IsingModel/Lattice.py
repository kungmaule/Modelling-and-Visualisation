import numpy as np
import sys
import random

##### DESCRIPTION ##############################################################

""" This file contains the code for the Lattice Object, which is used to run the
    Ising simulations. """

##### CLASS ####################################################################

class Lattice(object):

    def __init__(self, GridLength, Method, J, kT):

        """ This method initialises the lattice with given parameters, and
        generates an intiail random spin layout. """

        # Set Parameters
        self.GridLength = GridLength
        self.Method = Method
        self.J = J
        self.kT = kT

        # Initialise Random Spin Distribution
        self.Spins = np.zeros(shape=(GridLength, GridLength), dtype=float)
        for i in range(GridLength):
            for j in range(GridLength):
                r = random.random()
                if r < 0.5:
                    self.Spins[i,j] = -1.
                elif r >= 0.5:
                    self.Spins[i,j] = 1.

    def Glauber(self):

        """ This method performs a random spin flip, determines the energy change,
        and performs a Metropolis test to see if the spin should be flipped. """

        # Determine random spin to flip
        iTrial = np.random.randint(0,self.Spins.shape[0])
        jTrial = np.random.randint(0,self.Spins.shape[1])

        NewEnergy = 0.0
        #Apply pbc and sum energy terms
        NewEnergy += self.Spins[self.pbc(iTrial-1, 0),jTrial]
        NewEnergy += self.Spins[self.pbc(iTrial+1, 0),jTrial]
        NewEnergy += self.Spins[iTrial, self.pbc(jTrial-1, 1)]
        NewEnergy += self.Spins[iTrial, self.pbc(jTrial+1, 1)]
        DeltaE = NewEnergy*2*self.J*self.Spins[iTrial,jTrial]

        # Perform Metropolis Test if energy difference is negative
        if DeltaE <= 0:
            self.Spins[iTrial, jTrial] *= -1.0
        else:
            # Otherwise perform probability check for spin flip
            p = min(1, np.exp(-DeltaE/self.kT))
            r = np.random.random()
            if r <= p:
                self.Spins[iTrial,jTrial] *= -1.0

    def Kawasaki(self):

        """ This method performs swaps two different spins, determines the
        energy change, and performs a Metropolis test to see if the swap should happen. """

        # Select two random spins, ensuring that they are not the same
        i1, i2, j1, j2 = 0,0,0,0
        while i1 == i2 and j1==j2:
            i1 = np.random.randint(0,self.Spins.shape[0])
            i2 = np.random.randint(0,self.Spins.shape[0])
            j1 = np.random.randint(0,self.Spins.shape[1])
            j2 = np.random.randint(0,self.Spins.shape[1])

        # Pass if the same spin, otherwise work out energy change
        if self.Spins[i1,j1] == self.Spins[i2,j2]:
            pass
        else:
            DeltaE = 0.0
            # iterate over the two sites and work out the energy changes
            for Set in ((i1,j1),(i2,j2)):
                iTrial = Set[0]
                jTrial = Set[1]
                NewEnergy = 0.0
                NewEnergy += self.Spins[self.pbc(iTrial-1, 0),jTrial]
                NewEnergy += self.Spins[self.pbc(iTrial+1, 0),jTrial]
                NewEnergy += self.Spins[iTrial, self.pbc(jTrial-1, 1)]
                NewEnergy += self.Spins[iTrial, self.pbc(jTrial+1, 1)]
                DeltaE += NewEnergy*2*self.J*self.Spins[iTrial,jTrial]

            # Check if the spins are nearest neighbours and adjust accordingly
            if (i1==i2) and ((j1 == self.pbc(j2+1, 1)) or (j1 == self.pbc(j2-1, 1))):
                DeltaE -= 2.0*self.J
            elif (j1==j2) and ((i1 == self.pbc(i2+1, 0)) or (i1 == self.pbc(i2-1, 0))):
                DeltaE -= 2.0*self.J

            # Perform Metropolis Test if change in energy is negative
            if DeltaE <= 0.0:
                self.Spins[i1,j1] *= -1
                self.Spins[i2,j2] *= -1
            else:
                # Otherwise perform test based on probability metric
                p = min(1, np.exp(-DeltaE/self.kT))
                r = np.random.random()
                # Perform swap if tesy passes
                if r <= p:
                    self.Spins[i1,j1] *= -1
                    self.Spins[i2,j2] *= -1

    def Sweep(self):

        """ This method performs the sweep, iterating over all lattice points
        and running the relevant method for the chosen dynamics. """

        # Iterate over all grid points
        for n in range((self.Spins.shape[0]*self.Spins.shape[1])):

            # Check method and run appropriate dynamic method
            if self.Method == 'Glauber':
                self.Glauber()
            elif self.Method == 'Kawasaki':
                self.Kawasaki()

    def UpdateTemperature(self, T):
        """ This method updates the temperature of the Lattice. """
        self.kT = T

    def pbc(self, x, axis):
        """ This method utilises periodic boundary conditions for a specific axis. """
        return np.mod(x, self.Spins.shape[axis])

    def DetermineEnergy(self):
        """ This method determines the total energy of the lattice """
        SumTerm = 0
        for i in range(self.Spins.shape[0]):
            for j in range(self.Spins.shape[1]):
                Point = self.Spins[i,j]
                # i-1, j
                SumTerm += self.Spins[np.mod(i-1, self.Spins.shape[0]),j]*Point
                # i+1, j
                SumTerm += self.Spins[np.mod(i+1, self.Spins.shape[0]),j]*Point
                # i, j-1
                SumTerm += self.Spins[i,np.mod(j-1, self.Spins.shape[1])]*Point
                # i, j+1
                SumTerm += self.Spins[i,np.mod(j+1, self.Spins.shape[1])]*Point
        return 0.5*(-self.J * SumTerm)


    def DetermineMagnetisation(self):
        """ This method determines the absolute value of the total lattice magnetisation """
        return np.abs(np.sum(self.Spins))
