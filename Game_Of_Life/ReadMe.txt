Game of Life : Read Me Document

Please find in the following a description and explanation of my submission for this checkpoint.

USAGE: GameOfLife.py <Lattice Size> <Randomise Lattice (Y/N)> <Y to Animate/N not to Animate>

Some user error handling has been implemented to ensure correct inputs. The inputs are as follows:
- Inputting "Y" for randomise lattice will instantize the model with a random grid, and will run data collection for the stability histogram.
  Inputting "N" will provide a menu of options of oscillators and spaceships to choose from. Choosing an oscillator will simply display it.
  Choosing a spaceship will collect center of mass data for that spaceship. Any other inputs for the above will be requested according to the mode
  the simulation to be run in, and should be self explanatory.
- Picking a lattice size of N will generate an NxN lattice of cells. These will initially be randomly assigned as either alive (1) or dead (0).
- Inputting "Y" will cause the lattice to be animated and output to the user. Inputting "N" will skip this step.

Some values have been hard-coded.
- Stability Limit = 10
- Maximum Number of Sweeps = 5000
- Oscillation Limit = 100

A class was written for the Lattice, in order to efficiently update it. The code for the class is in Lattice.py.

A maximum limit of 5000 steps was given to prevent the simulation over-running (e.g. if it got caught into a high period oscillating state). For low period
oscillation states, an Oscillating attribute was used to keep track, and would abort the run if 100 steps of oscillation were detected (e.g. a beacon).
The stability limit chosen was ten steps. The simulation will complete when a user specified number of data points for stabilisation have been collected.
Histogram Data was recorded into a text file for ease of access (as it is a one dimensional list of values.)

For spaceships, the spaceships were run over a user specified number of steps and the center of mass recorded. CoM values were not taken if any cells were
adjacent to the boundary. The data is read back in in the FitVelocity.py file. A period of motion between two boundaries is found, and the velocity fit
for each dimension is found. If the spaceship moves diagonally, the speed is the speed in a single direction (e.g. Glider moves c/4 in x and y, so it's
overall speed is c/4.) The fitting was performed simply using a 1st order numpy polyfit.

Backup files of the all data files are saved in Data/Backup. These files were used to generate the graphs found in Plots. The plots were generated
by reading in the tables and plotting them using PlotHistogram.py, or for spaceship velocities / COM using FitVelocity.py.

The stabilisation Histogram can be found under plots.

The determined glider velocity was 0.25c in each direction. The overall magnitude of the velocity was 0.35c at 45 deg from the positive x axis.
The determined LWSS velocity was 0.5c in x, 0 in y. The overall magnitude of the velocity was 0.55c in the direction of the positive x axis.
The determined lobster velocity was -0.14c in each direction. The overall magnitude of the velocity was 0.2c at -135 deg from the positive x axis.
