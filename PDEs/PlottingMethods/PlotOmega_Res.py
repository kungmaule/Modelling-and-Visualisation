# Import relevant libraries
import numpy as np
import sys, random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from astropy.table import Table
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
rc('font', size=12)

##### MAIN #####################################################################

if __name__=="__main__":

    # Check inputs
    if len(sys.argv) == 1:
        print(f'Usage: python {sys.argv[0]} <Grid Size>')
    else:
        N = int(sys.argv[1])
        Field = sys.argv[2]

    n = 2000 if N==50 else 10000
    # Read in Data File
    SORData = Table.read(f'../Data/SOR/{Field}-Field/N{N}_ConvergenceData_Res.dat', format='ascii')
    Steps = SORData['steps']
    Omega = SORData['omega']
    arg = np.argmin(Steps)

    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(Omega[Steps<n], Steps[Steps<n], lw=.5, color='black')
    ax.plot(Omega[arg], Steps[arg], 'x', color='blue', label='Minimum $\omega$: %.4f' % Omega[arg])
    ax.set_xlabel('$\omega$')
    ax.set_ylabel('Steps')
    ax.set_xlim(1.95, 2.0)
    ax.legend(loc='upper left')
    ax.set_title('Plot of number of steps for convergence as a function of $\omega$')
    plt.savefig(f'../Plots/SOR/Convergence/N{N}_ConvergencePlot_Res_{Field}.png')
