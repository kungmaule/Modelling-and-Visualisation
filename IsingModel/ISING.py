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

""" This method generates the Energy, SHC, Magnetisation and susceptibility by
    iterating over a temperature range. """

##### METHODS ##################################################################

def StdError(Mean, MeanSqr, n):
    """ Determines the error on a quantity from its mean and sqr mean """
    return np.sqrt((MeanSqr - (Mean**2))/n)

def BootstrapError(Array, k, Func, Params):

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

        # Determine Mean and Mean Square
        Mean = np.mean(Bootstrap)
        MeanSqr = np.mean(Bootstrap*Bootstrap)

        # Determine the relevant value (c, x) and append it and it's square to arrays
        Val = Func(Mean, MeanSqr, *Params)
        V1 = np.append(V1, Val)
        V2 = np.append(V2, Val*Val)

    # Determine error via the averages of both arrays
    return np.sqrt(np.mean(V2) - (np.mean(V1)**2))

def Susceptibility(Mean, MeanSqr, N, kT):
    """ Determines the magnetic susceptibility """
    return (MeanSqr-(Mean**2))/(N*kT)

def SpecificHeatCapacity(Mean, MeanSqr, N, kT):
    """ Determines the specific heat capacity """
    return (MeanSqr - (Mean**2))/(N*(kT**2))

##### MAIN METHOD ##############################################################

if __name__=="__main__":

    # Check inputs for correct parameters
    if len(sys.argv) != 5:
        print(f'Incorrect Usage. \nCorrect Usage: python {sys.argv[0]} <G/K> <N> <Steps> <Animate(Y/N)>')
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
    if sys.argv[4].upper() == 'Y':
        Animate = True
    elif sys.argv[4].upper() == 'N':
        Animate = False
    else:
        print('Invalid Animate Input. \nPlease input "Y" to animate, or "N" for no animation.')
        quit()


    # Hard-coded Parameters
    J = 1.0
    InitialTemperature, FinalTemperature, Interval = 1.0, 3.0, 0.1
    N = GridLength * GridLength
    EqSteps = 100
    Tau = 10

    # Define the Lattice
    print('Initialising Lattice...')
    IsingLattice = Lattice(GridLength, Method, J, InitialTemperature)

    # Animate if Desired
    if Animate == True:
        fig = plt.figure(figsize=(6,6))
        col = cm.get_cmap('viridis', 4)
        newcol = ListedColormap(col(np.linspace(0.4,0.6,256)))
        im = plt.imshow(IsingLattice.Spins, animated=True, cmap=newcol)

    # Define Arrays for Temp, mean energy, mean magnetisation, shc, susceptibility
    Temperatures = np.arange(InitialTemperature, FinalTemperature+Interval, Interval)
    E_Arr = np.empty(Temperatures.shape[0])
    M_Arr = np.empty(Temperatures.shape[0])
    C_Arr = np.empty(Temperatures.shape[0])
    X_Arr = np.empty(Temperatures.shape[0])

    # Define arrays for errors
    E_Err = np.empty(Temperatures.shape[0])
    M_Err = np.empty(Temperatures.shape[0])
    C_Err = np.empty(Temperatures.shape[0])
    X_Err = np.empty(Temperatures.shape[0])

    # Iterate over all temperatures
    # Note here: use reversed, as we want to evolve the system from disordered to ordered
    for t in reversed(range(Temperatures.shape[0])):

        # Assign Temperature and update lattice
        T = Temperatures[t]
        print(f'Current Temperature: {T.round(2)}', end='\r')
        IsingLattice.UpdateTemperature(T)

        # Assign Temporary Arrays
        E = np.empty(0)
        M = np.empty(0)

        # Iterate over Number of Steps
        for n in range(StepNumber):

            # Perform Sweep
            IsingLattice.Sweep()

            # Allow 100 steps to equilibriate
            if n <= EqSteps:
                continue

            # Otherwise record data every 10th step
            if (n%Tau) == 0:

                # Animate if desired.
                if Animate == True:
                    plt.cla()
                    im=plt.imshow(IsingLattice.Spins, animated=True,cmap=newcol)
                    plt.draw()
                    plt.pause(0.0001)

                # Check Method
                if IsingLattice.Method != 'Kawasaki':

                    # Only Determine Magnetisation in Glauber Regime
                    Magnetisation = IsingLattice.DetermineMagnetisation()
                    M = np.append(M, Magnetisation)

                # Determine Lattice Energy
                Energy = IsingLattice.DetermineEnergy()
                E = np.append(E, Energy)

        # Calculate Means of Energy + Specific Heat
        MeanEnergy = np.mean(E)
        MeanSqrEnergy = np.mean(E*E)
        SHC = SpecificHeatCapacity(MeanEnergy, MeanSqrEnergy, N, T)

        # Update Arrays
        E_Arr[t] = MeanEnergy
        C_Arr[t] = SHC

        # Generate Error in E
        E_Err[t] = StdError(MeanEnergy, MeanSqrEnergy, E.shape[0])

        # Generate Error for Specific Heat Capacity and add to array
        C_Err[t] = BootstrapError(E, 1000, SpecificHeatCapacity, (N,T))

        # Write out to File
        EnergyTable = Table([Temperatures, E_Arr, E_Err, C_Arr, C_Err], names=['T', 'E', 'E_Err', 'SHC', 'SHC_Err'])
        EnergyTable.write(f'Data/EnergyData_{IsingLattice.Method}.dat', format='ascii', overwrite=True)

        # Check Method
        if IsingLattice.Method != 'Kawasaki':

            # Convert M into m by dividing by N
            #M = M/N

            # Calculate Means and susceptibility
            MeanMagnetisation = np.mean(M)
            MeanSqrMagnetisation = np.mean(M*M)
            Chi = Susceptibility(MeanMagnetisation, MeanSqrMagnetisation, N, T)

            # Update Arrays
            M_Arr[t] = MeanMagnetisation
            X_Arr[t] = Chi

            # Determine Error in M
            M_Err[t] = StdError(MeanMagnetisation, MeanSqrMagnetisation, M.shape[0])

            # Determine Error in Susceptibility
            X_Err[t] = BootstrapError(M, 1000, Susceptibility, (N,T))

            # Write out data to Table
            MagTable = Table([Temperatures, M_Arr, M_Err, X_Arr, X_Err], names=['T', 'M', 'M_Err', 'X', 'X_Err'])
            MagTable.write(f'Data/MagData_{IsingLattice.Method}.dat', format='ascii', overwrite=True)

    # Finish Simulation
    print('Simulation Completed.')
