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
        print(f'Usage: python {sys.argv[0]} <Method Name> <Grid Size> <Field>')
        quit()
    else:
        Root = sys.argv[1]
        GridSize = int(sys.argv[2])
        FieldID = sys.argv[3]
        Field = 'Electric' if FieldID == 'E' else 'Magnetic'
        color = '#00916E' if FieldID == 'E' else '#482677'

    # Set up Arrays
    X = np.arange(GridSize)
    Y = np.arange(GridSize)
    Phi = np.empty(shape=(GridSize, GridSize))
    Fx = np.empty(shape=(GridSize, GridSize))
    Fy = np.empty(shape=(GridSize, GridSize))
    Fz = np.empty(shape=(GridSize, GridSize))

    with open(f'../Data/{Root}/{FieldID}-Field/N{GridSize}_Slice.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split('       ')
            i = int(Tokens[0])
            j = int(Tokens[1])
            Phi[i,j] = float(Tokens[2])
            Fx[i,j] = float(Tokens[3])
            Fy[i,j] = float(Tokens[4])
            Fz[i,j] = float(Tokens[5])

    # Plot Phi
    fig = plt.figure(figsize=(8,8))
    plt.axis('off')
    im = plt.imshow(Phi, animated=False,interpolation='gaussian', cmap='magma')
    plt.colorbar(ticks=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25])
    plt.title(f'Plot of Plane {Field} Potential [{Root}]')
    plt.savefig(f'../Plots/{Root}/{FieldID}-Field/N{GridSize}_PotentialPlot.png')

    # Plot Field
    u, v = np.meshgrid(np.arange(-GridSize/2, GridSize/2), np.arange(-GridSize/2,GridSize/2))
    Normalisation = np.power(np.power(Fx,2) + np.power(Fy,2) + np.power(Fz,2), 1/2)
    Normalisation[Normalisation==0]=1
    fig, ax = plt.subplots(figsize=(10,10))
    ax.quiver(Y, X, Fz/Normalisation, Fy/Normalisation, color=color)
    ax.set_title(f'{Field} Field')
    plt.savefig(f'../Plots/{Root}/{FieldID}-Field/N{GridSize}_VectorPlot.png')
