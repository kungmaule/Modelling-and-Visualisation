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
        print(f'Usage: python {sys.argv[0]} <Grid Size>')
        quit()
    else:
        Root = sys.argv[1]
        GridSize = int(sys.argv[2])
        FieldID = 'E'
        Field = 'Electric'
        PotentialID = '\phi'

    # Set up Arrays
    X = np.arange(GridSize)
    Y = np.arange(GridSize)
    Phi = np.empty(shape=(GridSize, GridSize))
    Ex = np.empty(shape=(GridSize, GridSize))
    Ey = np.empty(shape=(GridSize, GridSize))
    Ez = np.empty(shape=(GridSize, GridSize))

    # Read in File
    with open(f'../Data/{Root}/{FieldID}-Field/N{GridSize}_Slice.txt', 'r') as file:
        for count, line in enumerate(file):
            Tokens = line.split('    ')
            i = int(Tokens[0])
            j = int(Tokens[1])
            Phi[i,j] = float(Tokens[2])
            Ex[i,j] = float(Tokens[3])
            Ey[i,j] = float(Tokens[4])
            Ez[i,j] = float(Tokens[5])

    R = np.sqrt((X-0.5*GridSize)**2 + (Y-0.5*GridSize)**2)
    E = np.zeros(R.shape[0])
    for i in range(len(E)):
        E[i] = np.sqrt(Ex[i,i]**2 + Ey[i,i]**2 + Ez[i,i]**2)
    for count, val in enumerate(R):
        if val == 0:
            break
        else:
            R[count] *= -1

    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(R[R!=0], E[R!=0], 'x', lw=.5, color='black')
    ax.axvline(0.0, linestyle='dashed', lw=.25, color='black')
    ax.set_xlabel('$r$')
    ax.set_ylabel('$%s(r)$' % FieldID)
    ax.set_title('Plot of $%s(r)$ [%s]' % (FieldID, Root))
    plt.savefig(f'../Plots/{Root}/{FieldID}-Field/N{GridSize}_Field.png')

    R_tot = np.array([np.sqrt((i-0.5*GridSize)**2+(j-0.5*GridSize)**2) for i in range(GridSize) for j in range(GridSize)])
    E_tot = np.array([np.sqrt(Ex[i,j]**2 + Ey[i,j]**2 + Ez[i,j]**2) for i in range(GridSize) for j in range(GridSize)])
    Phi_tot = np.array([Phi[i,j] for i in range(GridSize) for j in range(GridSize)])
    logR = np.log(R_tot[R_tot>0])
    logE = np.log(E_tot[R_tot>0])
    logPhi = np.log(Phi_tot[R_tot>0])

    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(logR, logE, 'x', color='black', label='Data Points')

    Mask = ((logR >= 1) & (logR <= 3))
    FitR = logR[Mask]
    FitE = logE[Mask]
    a, b = np.polyfit(FitR, FitE, 1)
    xfit = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    yfit = a*xfit + b
    ax.plot(xfit, yfit, color='red', lw=.5, label=f'Fit Gradient: {a:.3f}')
    ax.legend(loc='best')
    ax.set_xlabel('$\ln r$')
    ax.set_ylabel('$\ln %s(r)$' % FieldID)
    ax.set_title('')
    ax.set_title('Plot of $\ln %s(r)$ Fitting [%s]' % (FieldID, Root))
    plt.savefig(f'../Plots/{Root}/{FieldID}-Field/N{GridSize}_FieldFit.png')

    # Write out Data
    args = FitR.argsort()
    SortR = FitR[args]
    SortE = FitE[args]
    OutputTable = Table([SortR, SortE], names=['R', 'E'])
    OutputTable.write(f'../Data/{Root}/{FieldID}-Field/FitData/N{GridSize}_FieldFit.dat', format='ascii', overwrite=True)


    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(logR, logPhi, 'x', color='black', label='Data Points')

    Mask = ((logR >= 0) & (logR <= 2))
    FitR = logR[Mask]
    FitPhi = logPhi[Mask]
    m, c = np.polyfit(FitR, FitPhi, 1)
    xfit = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    yfit = m*xfit + c
    ax.plot(xfit, yfit, color='red', lw=.5, label=f'Fit Gradient: {m:.3f}')
    ax.legend(loc='best')
    ax.set_xlabel('$\ln r$')
    ax.set_ylabel('$\ln %s(r)$' % PotentialID)
    ax.set_title('')
    ax.set_title('Plot of $\ln %s(r)$ Fitting [%s]' % (PotentialID, Root))
    plt.savefig(f'../Plots/{Root}/{FieldID}-Field/N{GridSize}_PotentialFit.png')

    # Write out Data
    args = FitR.argsort()
    SortR = FitR[args]
    SortPhi = FitPhi[args]
    OutputTable = Table([SortR, SortPhi], names=['R', 'Phi'])
    OutputTable.write(f'../Data/{Root}/{FieldID}-Field/FitData/N{GridSize}_PotentialFit.dat', format='ascii', overwrite=True)
