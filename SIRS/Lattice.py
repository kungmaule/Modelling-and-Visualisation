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
    Ising simulations. """

##### CLASSES ####################################################################

class Lattice(object):

    """ Describe """

    ##### SET-UP METHODS #######################################################

    def __init__(self, GridLength, Animate=True, Immunity=False):

        """ This method initialises the lattice with given parameters, and
        generates an intiail random spin layout. """

        # S = 1 for Recovered
        # S = 0 for Susceptible
        # S = -1 for Infected
        # S = 2 for Immune

        # Set Parameters
        self.GridLength = GridLength
        self.N = self.GridLength * self.GridLength
        self.Animate = Animate
        self.Grid = np.zeros(shape=(GridLength, GridLength))
        self.Immunity = Immunity

        #self.SetInitialVariables()
        if self.Animate == True:
            self.InitaliseAnimation()


    def RandomiseLattice(self):
        # Determine fractions of states expected
        fInf, fSus, fRec = (1-self.fImm)/3, (1-self.fImm)/3, (1-self.fImm)/3
        # Randomise grid according to probabilities
        self.Grid = np.random.choice([-1, 0, 1, 2], size=(self.GridLength, self.GridLength), p=[fInf, fSus, fRec, self.fImm])

    def SetGrid(self, p1, p2, p3, fImm = 0):
        # Set all probabilities
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.fImm = fImm
        # Randomise the lattice according to the above probabilities
        self.RandomiseLattice()

    def InitaliseAnimation(self):
        fig = plt.figure(figsize=(6,6))
        plt.axis('off')
        # Check for Immunity before setting up colorbar
        if self.Immunity == False:
            self.cmap = ListedColormap(['#FDE725','#3CBB75', '#2D708E'])
            im = plt.imshow(self.Grid, animated=True, cmap=self.cmap, vmin=-1, vmax=1)
            plt.colorbar(ticks=[-1, 0, 1])
        else:
            self.cmap = ListedColormap(['#FDE725','#3CBB75', '#2D708E', '#481567'])
            im = plt.imshow(self.Grid, animated=True, cmap=self.cmap, vmin=-1, vmax=2)
            plt.colorbar(im, ticks=[-1, 0, 1, 2])

    def AnimateFrame(self):
        plt.cla()
        plt.axis('off')
        # Check for Immunity before setting up colourmap
        if self.Immunity == False:
            im=plt.imshow(self.Grid, animated=True,cmap=self.cmap, vmin=-1, vmax=1)
        else:
            im=plt.imshow(self.Grid, animated=True,cmap=self.cmap, vmin=-1, vmax=2)
        plt.draw()
        plt.pause(0.0001)

    ##### METHODS ##############################################################

    def Sweep(self):

        # Perform N
        for n in range(self.GridLength*self.GridLength):

            # Select a random cell + random value
            i,j = np.random.randint(0,self.Grid.shape[0]), np.random.randint(0,self.Grid.shape[1])
            r = np.random.random()

            # If state is susceptible, check for infected nearest neighbours and update accordingly
            if self.Grid[i,j] == 0:
                InfN = False
                for nn in [(self.pbc(i-1, 0), j), (self.pbc(i+1, 0), j), (i, self.pbc(j-1, 1)), (i, self.pbc(j+1, 1))]:
                    if self.Grid[nn] == -1:
                        InfN = True
                        break
                if InfN == True:
                    if self.p1 >= r:
                        self.Grid[i,j] -= 1

            # If state is recovered, update accordingly
            elif self.Grid[i,j] == 1:
                    if self.p3 >= r:
                        self.Grid[i,j] -= 1

            # If state is infected, update accordingly
            elif self.Grid[i,j] == -1:
                    if self.p2 >= r:
                        self.Grid[i,j] += 2


        # Update animation if required
        if self.Animate == True:
            self.AnimateFrame()

    def pbc(self, x, axis):
        """ This method utilises periodic boundary conditions for a specific axis. """
        return np.mod(x, self.Grid.shape[axis])

    def CountInfected(self):
        """ This method counts the number of infected cells """
        return np.count_nonzero(self.Grid == -1)

    def CountRecovered(self):
        """ This method counts the number of recovered cells """
        return np.count_nonzero(self.Grid == 1)

    def CountSusceptible(self):
        """ This method counts the number of susceptible cells """
        return np.count_nonzero(self.Grid == 0)

    def Variance(self, Infectivity):
        """ This method determines the variance in the infected fraction """
        Mean = np.mean(Infectivity)
        MeanSqr = np.mean(Infectivity**2)
        return (MeanSqr - Mean**2)/self.N

    def BootstrapError(self, Array, k=1000):

        """ Determines the error on a value via generating 1000 bootstrap samples,
        determining the quantity for each bootstrap array, generating means and Square
        means, from which an error is determined."""

        # Set Values for parameters
        V1 = np.empty(0)
        V2 = np.empty(0)

        # Iterate over k measurements
        for i in range(k):

            # Determine Bootstrap Array of same length as array
            Bootstrap = np.random.choice(Array, Array.shape[0], replace=True)

            # Determine the relevant value (c, x) and append it and it's square to arrays
            Val = self.Variance(Bootstrap)
            V1 = np.append(V1, Val)
            V2 = np.append(V2, Val*Val)

        # Determine error via the averages of both arrays
        return np.sqrt(np.mean(V2) - (np.mean(V1)**2))

    ##### DATA COLLECTION METHODS ##############################################

    def Run(self, StepNumber, StopIfDead):

        """ This method performs <StepNumber> runs of the simulation using it's
        current attributes as parameters. It allows for 100 equilibration steps
        It stops and returns [0] if the system reaches an absorbing state,
        otherwise it returns the mean infectivity fraction for the run. """

        # Define Constants
        EquilbrationSteps = 100

        # Define Arrays
        Infectivity = np.zeros(0)

        for step in range(StepNumber+EquilbrationSteps):

            # Sweep Lattice
            self.Sweep()

            # Count States
            Inf = self.CountInfected()
            Rec = self.CountRecovered()
            Sus = self.CountSusceptible()

            # Check site totals
            if Inf == 0 and StopIfDead == True:
                #print('\n System has reached an absorbing state.                    ')
                Infectivity = np.zeros(1)
                break

            # Allow Equilbration Steps
            if step <= EquilbrationSteps:
                #print(f'Step: {step}    Lattice Equilibrating', end='\r')
                pass

            else:

                # Check site totals
                if Inf == 0 and StopIfDead == True:
                    #print('\n System has reached an absorbing state.                    ')
                    Infectivity = np.zeros(1)
                    break

                else:
                    #print(f'Step: {step}            Taking Run Data', end='\r')
                    Infectivity = np.append(Infectivity, Inf)

        return Infectivity

    ##### Mode Codes ###########################################################

    def SingleRun(self, StepNumber, p1, p2, p3):

        """ This method runs a single run of the code, with animation, and outputs
        the infectivity fraction for the run against the steps. """

        # Mode
        Mode = 'Run' # Arbitrary

        # Run Lattice
        self.SetGrid(p1,p2,p3)

        # Run for Infectivity
        Infectivity = self.Run(StepNumber, StopIfDead = False)
        Infectivity /= self.N
        Steps = np.arange(1, StepNumber, 1)

        # Save Data
        OutputTable = Table([Steps, Infectivity], names=['Steps', 'Infectivity'])
        OutputTable.write(f'Data/{Mode}_Data.dat', format='ascii', overwrite=True)

    def Generate_p1p3_Plane(self, StepNumber, pInterval, p2=0.5):

        """ This method iterates over various p1, p3 values to generate data for
        a colourmap plane. It iterates from 0 to 1 in p1 and p3 in steps of pInterval
        recording the infectivity fraction and variance in an array, which is
        then saved and output. """

        # Hardcoded Variables
        pRange = np.arange(0.0, 1+pInterval, pInterval)
        ScaledAvgInf = np.zeros(shape=(pRange.shape[0], pRange.shape[0]))
        Variance = np.zeros(shape=(pRange.shape[0], pRange.shape[0]))

        for c1, p1 in enumerate(pRange):
            for c3, p3 in enumerate(pRange):

                # Generate Lattice
                self.SetGrid(p1,p2,p3)
                Infectivity = self.Run(StepNumber, StopIfDead=True)

                # Add data to arrays
                ScaledAvgInf[c1,c3] = np.mean(Infectivity)/self.N
                Variance[c1,c3] = self.Variance(Infectivity)

                print('\n', p1, p3, ' - Done')

        # Write out data
        names = [f'v{p}' for p in pRange]
        InfTable = Table(ScaledAvgInf, names=names)
        InfTable.write(f'Data/Plane_InfFrac[{pInterval}].dat', format='ascii', overwrite=True)
        VarTable = Table(Variance, names=names)
        VarTable.write(f'Data/Plane_Var[{pInterval}].dat', format='ascii', overwrite=True)

    def FindWaves(self, StepNumber, pInterval, p1i, p1f, p2, p3):

        """ This method fixes p2 and p3, and iterates over a fixed region of p1
        in steps, recording the variance and associated error. This data is
        then saved for the given conditions."""

        # Set up arrays
        pRange = np.arange(p1i, p1f+pInterval, pInterval)
        Variance = np.zeros(pRange.shape[0])
        VarianceError = np.zeros(pRange.shape[0])

        # Iterate over various p1
        for count, p1 in enumerate(pRange):

            # Set Lattice
            self.SetGrid(p1, p2, p3)
            Infectivity = self.Run(StepNumber, StopIfDead=True)

            # Record Data
            Variance[count] = self.Variance(Infectivity)
            VarianceError[count] = self.BootstrapError(Infectivity)

            print(f'[Completion: {round(100*count/len(pRange),2)} % ]', end='\r')

        # Save Data
        OutputTable = Table([pRange, Variance, VarianceError], names=['Probability', 'Var', 'Var_Err'])
        OutputTable.write(f'Data/Cut_p[{pInterval}].dat', format='ascii', overwrite=True)

    def ImmunityRun(self, StepNumber, fImmInterval, p1, p2, p3):

        """ This method runs an immune fraction run. Iterating from an immune
        fraction of 0 to 1 in steps of fImmInterval, it performs multiple runs
        per fImm in order to generate error bars, then records the mean and std
        error on the mean to .dat files. """

        # Assign Immune fraction range to check
        fImmRange = np.arange(0, 1+fImmInterval, fImmInterval)
        NumberOfRuns = 5

        # Assign arrays for data collection
        InfFrac = np.zeros(fImmRange.shape[0])
        Error = np.zeros(fImmRange.shape[0])

        # Iterate over the various immune fractions
        for count, fImm in enumerate(fImmRange):

            print(f'[Completion: {round(100*count/len(fImmRange),2)} % ]', end='\r')

            # Temporary array for the run repetitions
            Run_Infectivity = np.zeros(NumberOfRuns)

            for run in range(NumberOfRuns):

                self.SetGrid(p1, p2, p3, fImm=fImm)
                Infectivity = self.Run(StepNumber, StopIfDead=True)

                Run_Infectivity[run] = np.mean(Infectivity)/self.N

            # Take Mean of Infectivity
            InfFrac[count] = np.mean(Run_Infectivity)
            Error[count] = np.std(Run_Infectivity)/np.sqrt(Run_Infectivity.shape[0])

            with open(f'{N}_Output_Immune.txt', 'a') as file:
                file.write(f'{fImm} done \n')

        # Output Data
        OutputTable = Table([fImmRange, InfFrac, Error], names=['fImm', 'Inf', 'Err'])
        OutputTable.write(f'Data/ImmunityData_N{self.N**2}.dat', format='ascii', overwrite=True)
