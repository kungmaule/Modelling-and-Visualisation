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
    Cahn-Hilliard simulations in the CahnHilliard class, as well as the
    code for running all variations of the Poisson/Magnetic simulations
    using boundary values. """

##### CLASSES ####################################################################

class CahnHilliard(object):

    ##### SET-UP METHODS #######################################################

    def __init__(self, GridLength, Phi0, Noise, a, M, k, Animate=False):

        """ This method initialises the lattice with given parameters, and
        generates an intiail random layout. """

        # Set simulation parameters
        self.GridLength = GridLength
        self.Animate = Animate
        self.deltat, self.deltax = 1, 1
        self.cmap = 'plasma'
        # Set initial parameters and generate lattice.
        self.Set(Phi0, Noise, a, M, k)

        # Initialise Animation if desired
        if self.Animate == True:
            self.InitaliseAnimation()

        # For recording screenshots of behaviour
        self.Record = False

        # For recording data
        self.Label = 'PhiHalf' if abs(Phi0) == 0.5 else 'PhiZero'

    def InitialiseLattice(self):

        """ This method generates a lattice with normally distributed noise
            based on the desired noise in Phi. """

        self.Phi = np.random.normal(loc=self.Phi0, scale=self.Noise, size=(self.GridLength, self.GridLength))


    def Set(self, Phi0, Noise, a, M, k):
        """ Sets lattice with running Parameters """
        self.Phi0 = Phi0
        self.Noise = Noise
        self.a = a
        self.M = M
        self.k = k
        # Initialises lattice
        self.InitialiseLattice()

    def InitaliseAnimation(self):
        """ This method initialises the animation of a single plane. """
        fig = plt.figure(figsize=(6,6))
        plt.axis('off')
        im = plt.imshow(self.Phi, animated=True,interpolation='gaussian', cmap=self.cmap)
        plt.colorbar()

    def AnimateFrame(self):
        """ This method animates an individual plane frame of the potential. """
        plt.cla()
        plt.axis('off')
        im = plt.imshow(self.Phi, animated=True, interpolation='gaussian', cmap=self.cmap)
        plt.draw()
        plt.pause(0.0001)

    def RecordPlane(self, step):

        """ For a given step, generates a mod value for recording the plane, and
            if the mod of the step is zero, save the plane."""

        if 0 <= step <= 10000:
            Mod = 1000
        elif 10000 < step <= 100000:
            Mod = 10000
        elif 100000 < step < 1000000:
            Mod = 100000
        else:
            Mod = 1000000
        if step % Mod == 0 and self.Save == True:
            # Record plane
            fig = plt.figure(figsize=(6,6))
            plt.axis('off')
            im = plt.imshow(self.Phi, animated=False, interpolation='gaussian', cmap=self.cmap, vmin=-1, vmax=1)
            plt.title('Cahn Hilliard: $\phi_{0} = %.1f$ [Step %s] ' % (self.Phi0, str(step)))
            plt.colorbar()
            plt.savefig(f'Plots/Cahn-Hilliard/SavedPlots/{self.Label}/Plane_{ID}.png')
            plt.close()
            # Increment planeID for next save
            self.PlaneID += 1


    ##### METHODS ##############################################################

    def RollArray(self, Array):

        """ This method rolls the array in 2D, returning a roll array in
            each dimension. """

        return np.roll(Array, -1, axis=0), np.roll(Array, 1, axis=0), np.roll(Array, -1, axis=1), np.roll(Array, 1, axis=1)

    def DiscreteLaplacian(self, Array):

        """ This method calculates the laplacian of the 2D array via rolling. """

        Arr_Up, Arr_Down, Arr_Right, Arr_Left = self.RollArray(Array)
        return (Arr_Up + Arr_Down + Arr_Left + Arr_Right - 4*Array)

    def Sweep(self):

        """ This method performs an individual sweep update. Firstly, mu is
            calculated, and then Phi is updated using the updated mu. """

        # Calculate Mu
        self.Mu = -(self.a*self.Phi) + (self.a*np.power(self.Phi, 3)) - ((self.k*self.DiscreteLaplacian(self.Phi))/(self.deltax**2))
        # Update Phi
        self.Phi += (self.M * self.deltat * self.DiscreteLaplacian(self.Mu))/(self.deltax**2)
        # Update animation if required
        if self.Animate == True:
            self.AnimateFrame()

    def Determinef(self):

        """ This method determines the free energy of the lattice at a given phi. """

        # Take gradient
        Grad_i, Grad_j = np.gradient(self.Phi)
        # Return free energy
        return (self.deltat**2) *  np.sum(-(self.a * 0.5 * self.Phi**2) + (self.a * 0.25 * self.Phi **4) + (self.k * 0.5 * (Grad_i**2+Grad_j**2)))

    ##### DATA COLLECTION METHODS ##############################################

    def Run(self, StepNumber):

        """ This method performs an actual simulation run, recording the free
            energy of the lattice when iterating over a number of steps. if
            record is true, will also take screenshots of Phi at various
            timesteps. FreeEnergyDensity and steps are returned to the main. """

        # Set up array to store data
        FreeEnergyDensity = np.zeros(0)
        Steps = np.zeros(0)

        for step in range(StepNumber):

            # Record free energy density every ten steps
            if step % 10 == 0:
                FreeEnergyDensity = np.append(FreeEnergyDensity, self.Determinef())
                Steps = np.append(Steps, step)
                print(f'[ Step: {step}/{StepNumber} - Completion: {round(100*step/StepNumber,3)} % ]', end='\r')

            # Sweep Lattice
            self.Sweep()

        return FreeEnergyDensity, Steps

    def TestParameters(self, Interval):

        """ This method performs a dt dx parameter test, looking for convergence.
            It produces an array of values. Steps to divergence (when any phi
            array element becomes nan) are recorded. If the simulation converges
            within 200 steps, a nan is instead recorded. """

        # Steps to convergence
        MaxSteps = 200
        # ID for ordering planes
        self.PlaneID = 0
        # Generate dx, dt array
        dx, dt = np.arange(Interval, 1+Interval, Interval), np.arange(Interval, 1+Interval, Interval)
        # Arrays for recording data
        steps, xstep, tstep = np.empty(0), np.empty(0), np.empty(0)

        # Iterate over all combinations
        for i, x in enumerate(dx):
            for j, t in enumerate(dt):
                # Update lattice parameters
                self.deltat = t
                self.deltax = x
                # Reset Lattice
                self.InitialiseLattice()
                # Add x, t values
                xstep = np.append(xstep, x)
                tstep = np.append(tstep, t)

                # Iterate over maximum sweeps
                for n in range(MaxSteps):
                    self.Sweep()
                    # If any lattice element becomes nan, record the steps and leave loop
                    if np.isnan(self.Phi).any():
                        steps = np.append(steps, n)
                        break
                # If max sweeps reached, add nan to steps.
                if n == MaxSteps-1:
                    steps = np.append(steps, np.nan)
                print(f'dx: {x:.2f}/{dx[-1]:.2f}    dt: {t:.2f}/{dt[-1]:.2f}', end='\r')

        return xstep, tstep, steps

class Poisson(object):

    ##### SET-UP METHODS #######################################################
    def __init__(self, GridLength, Phi0, Noise, i, j, k, Field, Animate=False):

        """ This method initialises the lattice with given parameters, and
            generates an intiail random potential layout. It also initialises
            the animation if desired.  """

        # Set out lattice parameters and simulation parameters
        self.GridLength = GridLength
        self.Animate = Animate
        self.deltat, self.deltax = 1, 1
        self.Field = Field
        """ IF MAGNETIC CHANGE BOUNDARY CONDITIONS TO USE PBC AT THE TOP AND BOTTOM OR WORK IN 2D"""
        # Generate lattice with given parameters
        self.Set(Phi0, Noise, i, j, k)
        # Set plane to be the plane of charge
        self.Plane = k
        # Set convergence limit
        self.Limit = 1e-3
        # Set boolean for convergence
        self.Stop = False
        # Initialise animation
        if self.Animate == True:
            self.InitaliseAnimation()


    def InitialiseLattice(self):
        """ This method initialises a lattice based on random noise, and
            enforces the Dirichlet BC. """
        # Initialise lattice with Phi0
        self.Phi = np.full(shape=(self.GridLength, self.GridLength, self.GridLength), fill_value=self.Phi0) + self.Noise*np.random.rand(self.GridLength, self.GridLength, self.GridLength)
        # Set boundaries to zero (Dirichlet)
        self.DirichletBC(self.Phi)

    def AddRho(self,i,j,k):
        """ This method either adds a monopole charge or a wire in the z
            direction based on given coordinates. """
        # Initialise empty lattice
        self.Rho = np.zeros(shape=(self.GridLength, self.GridLength, self.GridLength))
        # Add Charge / Wire
        if self.Field == 'E':
            self.Rho[i,j,k] = 1
        else:
            self.Rho[i,j,:] = 1

    def Set(self, Phi0, Noise, i, j, k):
        """ Sets lattice with running Parameters, initialises and adds Rho. """
        self.Phi0 = Phi0
        self.Noise = Noise
        self.InitialiseLattice()
        self.AddRho(i,j,k)

    def DirichletBC(self, Array):
        """ This method enforces Dirichlet boundary conditions, setting the
            ends of the simulation box to zero whenever called. """
        for Boundary in [0, -1]:
            Array[Boundary,:,:] = 0
            Array[:,Boundary,:] = 0
            Array[:,:,Boundary] = 0

    def InitaliseAnimation(self):
        """ This method initialises the animation for a given plane potential. """
        fig = plt.figure(figsize=(6,6))
        plt.axis('off')
        self.cmap = 'magma'
        im = plt.imshow(self.Phi[:,:,self.Plane], animated=True,interpolation='gaussian', cmap=self.cmap)
        plt.colorbar()

    def AnimateFrame(self):
        """ This method animates an individual plane frame of the potential. """
        plt.cla()
        plt.axis('off')
        im = plt.imshow(self.Phi[:,:,self.Plane], animated=True, interpolation='gaussian', cmap=self.cmap)
        plt.draw()
        plt.pause(0.0001)

    ##### METHODS ##############################################################

    def GenerateMask(self):

        """ This method generates a checkerboard mask array of 1s and 0s
            for use in the Gauss Seidel algorithm. """

        # Make array of zeros
        Mask = np.zeros(shape=(self.GridLength, self.GridLength, self.GridLength))
        # Set alternating elements to 1
        Mask[::2, ::2, 1::2] = 1
        Mask[::2, 1::2, ::2] = 1
        Mask[1::2, 1::2, 1::2] = 1
        Mask[1::2, ::2, ::2] = 1
        # Store as object property
        self.Mask = Mask

    def RollArray3D(self, Array):

        """ This method performs a 3D roll of a given array, returning a rolled
            array in each direction. """

        return np.roll(Array, -1, axis=0), np.roll(Array, 1, axis=0), np.roll(Array, -1, axis=1), np.roll(Array, 1, axis=1), np.roll(Array, -1, axis=2), np.roll(Array, 1, axis=2)

    def ElectricField(self):

        """ Determine the electric field from the - of the gradient of
            the potential """

        return -1 * np.array(np.gradient(self.Phi))

    def MagneticField(self):

        """ Determine the magnetic field assuming Ax, Ay are zero. """

        # Take gradients of the field and return them in the correct order
        Gradients = np.array(np.gradient(self.Phi))
        return np.array([Gradients[1], -Gradients[0], Gradients[2]])

    def Jacobi(self):

        """ This method performs a single Jacobi sweep, rolling the phi array
            and updating phi based on the neighbours. BC are enforced. """

        Arr_Up, Arr_Down, Arr_Right, Arr_Left, Arr_Forward, Arr_Back = self.RollArray3D(self.Phi)
        self.Phi = (Arr_Up + Arr_Down + Arr_Left + Arr_Right + Arr_Forward + Arr_Back + np.power(self.deltax,2)*self.Rho)/6.
        # Enforce PBC
        self.DirichletBC(self.Phi)

    def GaussSeidel(self, Array, Mask):

        """ This method performs a Gauss Seidel update. For b and w squares
            on the given mask, the lattice is rolled and the respective site
            updated, with boundary conditions enforced.

            This method has a return for use in the SOR algorithm   """

        # Iterate over both sublattices
        for bw in [0, 1]:
            # Roll Arrays
            Arr_Up, Arr_Down, Arr_Left, Arr_Right, Arr_Forward, Arr_Back = self.RollArray3D(Array)
            # Update relevant sublattice
            Array[Mask==bw] = (1/6)*(Arr_Up[Mask==bw] + Arr_Down[Mask==bw] + Arr_Left[Mask==bw] + Arr_Right[Mask==bw] + Arr_Forward[Mask==bw] + Arr_Back[Mask==bw] + (self.deltax**2)*self.Rho[Mask==bw])
            # Reinforce Boundary Conditions
            self.DirichletBC(Array)
        return Array

    def SOR(self, Mask, omega):

        """ This method performs a single sweep of SOR< using omega and G-S.
            All enforcement of boundary conditions is handled by self.GaussSeidel """
        """ THIS NEEDS TO WORK ON THE SUBLATTICES I THINK? HENCE NOT PROPER OVER RELAX EXPECT OMEGA = 1.94"""
        return (1-omega)*self.Phi + omega*self.GaussSeidel(self.Phi, Mask)

    def WriteOutSliceData(self, Method):

        """ This method writes out the data for the potential and Field
            on a single slice given by self.Plane to data files. """

        # Get Field Data in specific plane
        xSlice = self.xField[:,:,self.Plane]
        ySlice = self.yField[:,:,self.Plane]
        zSlice = self.zField[:,:,self.Plane]

        # Write out Plane Data
        with open(f'Data/{Method}/{self.Field}-Field/N{self.GridLength}_Slice.txt', 'a') as file:
            for i in range(self.GridLength):
                for j in range(self.GridLength):
                    file.write(f'{i:3}    {j:3}    {self.Phi[i,j,self.Plane]:.12f}    {xSlice[i,j]:.12f}    {ySlice[i,j]:.12f}    {zSlice[i,j]:.12f}\n')

        # Write out Full Data Set (only if not too big)
        if self.GridLength <= 50:
            with open(f'Data/{Method}/{self.Field}-Field/N{self.GridLength}_Full.txt', 'a') as file:
                for i in range(self.GridLength):
                    for j in range(self.GridLength):
                        for k in range(self.GridLength):
                            file.write(f'{i:3}    {j:3}    {k:3}    {self.Phi[i,j,k]:.8f}    {self.xField[i,j,k]:.8f}    {self.yField[i,j,k]:.8f}    {self.zField[i,j,k]:.8f}\n')

    ##### DATA COLLECTION METHODS ##############################################

    def JacobiRun(self):

        """ This method runs a single iteration of the jacobi algorithm
            until convergence, and records field and potential data. """

        # Step Variable
        step = 0

        while self.Stop == False:

            # Store old lattice
            self.PhiOld = np.copy(self.Phi)

            # Sweep
            self.Jacobi()

            # Check for Convergence
            Difference = self.Phi - self.PhiOld
            self.PhiDifference = np.sum(np.abs(Difference))
            self.Stop = True if self.PhiDifference <= self.Limit else False

            # Animate
            if self.Animate == True:
                self.AnimateFrame()
            step += 1
            # Output
            print(f'Step: {step}    Δφ: {self.Limit} < {self.PhiDifference:.6f}', end='\r')

        print('     -=-=-=-=- Simulation Completed -=-=-=-=-    ')

        # Determine the relevant field
        Field = self.ElectricField() if self.Field == 'E' else self.MagneticField()
        self.xField, self.yField, self.zField = Field[0], Field[1], Field[2]

        # Record Data
        self.WriteOutSliceData('Jacobi')

    def GaussSeidelRun(self):

        """ This method runs a single iteration of the Gauss Seidel
            algorithm until convergence, and records field and potential data.    """

        # Generate Mask
        self.GenerateMask()
        step = 0

        while self.Stop == False:

            # Copy old Array
            self.PhiOld = np.copy(self.Phi)

            # Sweep Lattice
            self.Phi = self.GaussSeidel(self.Phi, self.Mask)

            # Determine Difference and stop if at limit
            Difference = self.Phi - self.PhiOld
            self.PhiDifference = np.sum(np.abs(Difference))
            self.Stop = True if self.PhiDifference <= self.Limit else False

            # Output to terminal + animate
            step += 1
            if self.Animate == True:
                self.AnimateFrame()
            print(f'Step: {step}    Δφ: {self.Limit} < {self.PhiDifference:.6f}', end='\r')

        print('     -=-=-=-=- Simulation Completed -=-=-=-=-    ')

        # Determine the Field
        Field = self.ElectricField() if self.Field == 'E' else self.MagneticField()
        self.xField, self.yField, self.zField = Field[0], Field[1], Field[2]
        # Write Out Data
        self.WriteOutSliceData('GaussSeidel')

    def SORRun(self, omega):

        """ This method runs a single iteration of the SOR algorithm until
            convergence, and records field and potential data. """

        # Iterate until convergence
        step = 0
        while self.Stop == False:

            # Copy old Array
            self.PhiOld = np.copy(self.Phi)

            # Sweep Lattice
            self.Phi = self.SOR(self.Mask, omega)

            # Determine Difference and stop if at limit
            Difference = self.Phi - self.PhiOld
            self.PhiDifference = np.sum(np.abs(Difference))
            self.Stop = True if self.PhiDifference <= self.Limit else False

            # Output to terminal + animate
            step += 1
            if self.Animate == True:
                self.AnimateFrame()
            print(f'Step: {step}    omega: {omega:.2f}      Δφ: {self.Limit} < {self.PhiDifference:.6f}', end='\r')
            # Handle case of divergence
            if step == self.MaxSteps:
                print('Failure to Converge')
                return step
        print('     -=-=-=-=- Simulation Completed -=-=-=-=-    ')
        return step

    def SOR_Single(self, omega):

        """ This method performs a single run using SOR at a given omega value,
            and records the field data """
        # Generate Mask
        self.GenerateMask()
        self.MaxSteps = 100000

        # Run - resultant steps are not necessary
        _ = self.SORRun(omega)

        # Record Values
        Field = self.ElectricField() if self.Field == 'E' else self.MagneticField()
        self.xField, self.yField, self.zField = Field[0], Field[1], Field[2]
        # Write Out Data
        self.WriteOutSliceData('SOR')

    def SOR_Convergence(self, Interval = 0.01):

        """ This method examines the convergence of the SOR algorithm for a
            varying omega with a given interval."""

        # Generate Mask
        self.GenerateMask()
        # Get number of values
        vals = int(1/Interval)
        # Study restricted range?
        Restricted = False
        # Generate Range of omega (including restricted more thorough range)
        Omega_Range = np.linspace(1, 2, vals, endpoint=False) if Restricted == False else np.linspace(1.9, 2.0, vals, endpoint=False)
        Label = '_Res' if Restricted == True else ''
        Steps = np.empty(Omega_Range.shape[0])
        self.MaxSteps = 20000 # hardlimit in case of divergence

        # iterate over omega range
        for count, omega in enumerate(Omega_Range):

            # Run at given omega
            Steps[count] = self.SORRun(omega)

            # Reset Lattice
            self.InitialiseLattice()
            self.Stop = False
            #print('.                                                    ', end='\r')

        # Record Data in Astropy Table
        SORTable = Table([Omega_Range, Steps], names=['omega', 'steps'])
        SORTable.write(f'Data/SOR/{self.Field}-Field/N{self.GridLength}_ConvergenceData{Label}.dat', format='ascii', overwrite=True)
