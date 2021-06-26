Ising Model: Read Me Document

Please find in the following a description and explanation of my submission for this checkpoint.

USAGE: ISING.py <G for Glauber/K for Kawasaki> <Lattice Size> <Number of Steps> <Y to Animate/N not to Animate>

Some user error handling has been implemented to ensure correct inputs. The inputs are as follows:
- Inputting "G" will instantize the model to use Glauber dynamics. Inputting "K" will similarly mean the model uses Kawasaki dynamics.
- Picking a lattice size of N will generate an NxN lattice of spins. These will initially be randomly assigned as either spin up or spin down.
- The number of steps assumed for all graphs/data outputs is 10000, though this has been left up to user input.
- Inputting "Y" will cause the lattice to be animated and output to the user. Inputting "N" will skip this step.

Some values have been hard-coded.
- 100 steps has been assumed to be the equilibration time.
- 10 steps has been assumed to be the autocorrelation time.
- The Lattice is assumed to be square, such that the number of particles is NxN.

A class was written for the Lattice, in order to efficiently update it. The code for the class is in lattice.py.

Note that the loop over different temperatures runs backwards. This was done to start the system in a more disordered state, before progressing
to an ordered state. This way, it is more likely that the system will avoid any meta-stable states that can form at low temperatures if initialised from
random. These states will be seen to be horizontal/vertical/diagonal lines of states in one spin configuration.

Data was output using astropy's tables - chosen as they are an efficient method of writing and reading arrays from .dat files. ISING.py will
produce a datafile for energy and specific heat capacity for each method, stored as 'EnergyData_(Method).dat', which has columns for the temperature
the reading was taken at, the mean energy, error on the mean energy, specific heat capacity, and error on specific heat capacity. The error in Energy
is determined via a standard error on the mean, the error in heat capacity is formed by 1000 bootstrap samples. 

If using Glauber dynamics there will also be a file generated for magnetisation data: 'MagData_Glauber.dat' which stores temperature, magnetisation,
error on magnetisation, susceptibility and error on susceptibility. The errors are calculated in the same manner as described for energy. These values
are not calculated for Kawasaki dynamics as the overall magnetisation of the lattice remains constant, because spins swap and do not individually flip.

Backup files of the three data files are saved in Data/Backup. These files were used to generate the graphs found in Plots. The plots were generated
by reading in the tables and plotting them using PlotFigures.py.

I have also included a version of the simulation to run simply to see the animation: ISING_Animate.py.
USAGE: ISING_Animate.py <G/K> <Lattice Size> <Number of Steps> <Temperature>
This version of the method simply has the loop for iterating over the temperatures removed and all data writing removed. It's simply to be able to
view the lattice at a specific temperature, so will not effect data files.

For both methods, the determined critical temperature was T=2.3, chosen as the maximum of the SHC or Susceptibility. The critical temperature has been
plotted on all the plots as a faint green line.
