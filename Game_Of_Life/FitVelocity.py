# Import Relevant Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys, time
from Lattice import Lattice
from astropy.table import Table
from matplotlib.ticker import PercentFormatter
from matplotlib import rc
rc('text', usetex=True)
from scipy.signal import argrelmax, argrelmin

##### DESCRIPTION ##############################################################



##### MAIN METHOD ##############################################################

if __name__=="__main__":

    print('Select a State to plot. Available States: ')
    PreCodedStates = ['Glider', 'Lwss', 'Lobster']
    for i, n in enumerate(PreCodedStates):
        print(f'{i}: {n}')

    i = None

    while i not in np.arange(0, len(PreCodedStates)):
        i = int(input('Choice: '))
    State = PreCodedStates[i]

    # Read in Data
    Data = Table.read(f'Data/CenterOfMass({State}).dat', format='ascii')
    Time = Data['Time']
    X = Data['X']
    Y = Data['Y']
    R = np.sqrt(X**2+Y**2)

    # Plot Position Graph
    plt.figure(figsize=(8,8))
    plt.plot(Time, R, color='black')
    plt.title(f'Plot of {State} position against steps.')
    plt.xlabel('Steps')
    plt.ylabel('COM Distance from Origin')
    plt.savefig(f'Plots/Position/{State}_Position_vs_Time.png')

    # Find Minimum and Maximum
    MinR = np.min(R)
    MaxR = np.max(R)
    MinN = np.where(R==MinR)[0][0]
    MaxN = np.where(R==MaxR)[0][1]

    if State == 'Lobster':
        MinN = 269 # Arbitrarily Chosen as Lobster moves weird
        MaxN = 367 # Arbitrarily Chosen as Lobster moves weird

    ContTimes = Time[MinN:MaxN]
    ContX = X[MinN:MaxN]
    ContY = Y[MinN:MaxN]

    plt.figure()
    plt.plot(ContTimes, ContX, 'x', label='X', color='red')
    plt.plot(ContTimes, ContY, 'x', label='Y', color='blue')
    ax, bx = np.polyfit(ContTimes, ContX, 1)
    print(f'X-Velocity: {ax.round(2)}c')
    ay, by = np.polyfit(ContTimes, ContY, 1)
    print(f'Y-Velocity: {ay.round(2)}c')
    plt.plot(ContTimes, ax*ContTimes + bx, color='red', alpha=0.5, label=f'X-Fit with v={ax.round(2)}')
    plt.plot(ContTimes, ay*ContTimes + by, color='blue', alpha=0.5, label=f'Y-Fit with v={ay.round(2)}')
    plt.xlabel('Steps')
    plt.ylabel('Position [X/Y]')
    plt.title(f'Plot of Velocity Fit for {State}')
    plt.legend(loc='best')
    plt.savefig(f'Plots/Velocity/{State}_VelocityFit.png')

    # Output Overall Velocity
    Magnitude = np.sqrt(ax**2 + ay**2)
    Direction = np.degrees(np.arctan2(ay,ax))

    print(f'Overall Velocity Magnitude {Magnitude.round(2)} in Direction {Direction.round(0)} degrees from x-axis.')
