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

    if len(sys.argv) == 1:
        print(f'Usage: python {sys.argv[0]} <Phi0> <Interval>')
        quit()
    else:
        Phi0 = float(sys.argv[1])
        Interval = float(sys.argv[2])

    ParamData = Table.read(f'../Data/Cahn-Hilliard/ParameterData/Phi0[{Phi0}].dat', format='ascii')
    dx = ParamData['dx']
    dt = ParamData['dt']
    steps = ParamData['s']
    Range = np.arange(Interval, 1+Interval, Interval)
    Matrix = np.array([steps[dx==i] for i in Range])
    # Plot Phi
    fig = plt.figure(figsize=(8,8))
    im = plt.imshow(Matrix, cmap='viridis', origin='lower', extent=[0,1,0,1],)
    cb = plt.colorbar()
    cb.set_label('Steps to Divergence (white for convergence)')
    plt.xlabel('dx')
    plt.ylabel('dt')
    plt.title(f'Plot of Parameter convergence')
    plt.savefig(f'../Plots/Cahn-Hilliard/Parameters/Phi0[{Phi0}]_Params.png')
    plt.show()
