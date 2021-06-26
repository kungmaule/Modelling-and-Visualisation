import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
rc('text', usetex=True)

##### MAIN #####################################################################

if __name__=="__main__":

    if len(sys.argv) == 1:

        # Determine which Figure to PLot
        PlotModes = ['Single Run', 'p1p3 Plane', 'p1 Cut', 'Immmunity Fraction']
        print('Select a plot to produce: ')
        for n, m in enumerate(PlotModes):
            print(f'{n+1}: {m}')
        Mode = -1
        while Mode not in np.arange(1, len(PlotModes)+1, 1):
            Mode = int(input('Number to plot: '))

    else:
        Mode = int(sys.argv[1])

    print('\n')

    if Mode == 1: # SINGLE RUN

        # Get right mode to plot
        PlotRuns = ['Dynamic Equilibrium', 'Absorbing State', 'Waves', 'All']
        print('Select a behaviour to plot: ')
        for n, m in enumerate(PlotRuns):
            print(f'{n+1}: {m}')
        Run = -1
        while Run not in np.arange(1, len(PlotRuns)+1, 1):
            Run = int(input('Number to plot: '))

        # Plot individual modes
        if 1 <= Run <= 3:
            if Run == 1:
                RunData = Table.read('Data/Dynamic_Equilibrium_Data.dat', format='ascii')
                Label = 'DE'
            elif Run == 2:
                RunData = Table.read('Data/Absorbing_State_Data.dat', format='ascii')
                Label = 'AS'
            elif Run == 3:
                RunData = Table.read('Data/Waves_Data.dat', format='ascii')
                Label = 'WA'

            Steps = RunData['Steps']
            Infectivity = RunData['Infectivity']

            plt.figure(figsize=(12, 6))
            plt.title(f'Infectivity vs Steps: {PlotRuns[Run-1]} Behaviour')
            plt.plot(Steps, Infectivity, color='black', lw=.25)
            plt.xlabel('Steps')
            plt.xlim(Steps[0], Steps[-1])
            plt.ylabel('Infectivity')
            plt.tight_layout()
            plt.savefig(f'Plots/Behaviour/{Label}_Plot.png')

        # Plot all modes overlayed
        else:
            DEData = Table.read('Data/Dynamic_Equilibrium_Data.dat', format='ascii')
            ASData = Table.read('Data/Absorbing_State_Data.dat', format='ascii')
            WAData = Table.read('Data/Waves_Data.dat', format='ascii')

            plt.figure(figsize=(12, 6))
            plt.title(f'Infectivity vs Steps: All Behaviours')
            plt.plot(DEData['Steps'], DEData['Infectivity'], color='blue', lw=.25, label='Dynamic Equilibrium')
            plt.plot(ASData['Steps'], ASData['Infectivity'], color='red', lw=.25, label='Absorbing State')
            plt.plot(WAData['Steps'], WAData['Infectivity'], color='green', lw=.25, label='Waves')
            plt.xlabel('Steps')
            plt.xlim(0, 5000)
            plt.ylabel('Infectivity')
            #plt.tight_layout()
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1), ncol=3)
            plt.savefig(f'Plots/Behaviour/All_Plot.png')

    elif Mode == 2: # PLANE CODE

        # Get Correct Interval
        Intervals = [0.2, 0.01, 0.05]
        pInterval = None
        while pInterval not in Intervals:
            print(f'Available Intervals: {Intervals[0]}, {Intervals[1]}, {Intervals[2]}')
            pInterval = float(input('Chosen Interval: '))

        # Set up Arrays
        pRange = np.arange(0, 1+pInterval, pInterval)
        FracInfArr = np.empty(shape=(0, pRange.shape[0]))
        VarArr = np.empty(shape=(0, pRange.shape[0]))

        # Read in Data
        FracInfData = Table.read(f'Data/Plane_InfFrac[{pInterval}].dat', format='ascii')
        VarData = Table.read(f'Data/Plane_Var[{pInterval}].dat', format='ascii')

        # Assign to Array
        for pVal in pRange:
            FracInfArr = np.vstack((FracInfArr, FracInfData[f'v{pVal}']))
            VarArr = np.vstack((VarArr, VarData[f'v{pVal}']))

        # Plot Fractional Infectivity Plane
        plt.figure(figsize=(8,8))
        plt.imshow(FracInfArr, cmap='viridis', origin='lower', extent=[0,1,0,1], vmin=0, vmax=0.5)
        plt.xlabel('Probability $p_{1}$')
        plt.ylabel('Probability $p_{3}$')
        plt.title('Heatmap of Fractional Infectivity for $p_{2} = 0.5$')
        plt.colorbar(ticks=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5],fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.savefig(f'Plots/Planes/Plane_InfFrac[{pInterval}].png')

        # Plot Variance Plane
        plt.figure(figsize=(8,8))
        plt.imshow(VarArr, cmap='viridis', origin='lower', extent=[0,1,0,1], vmin=0, vmax = 1.4)
        plt.xlabel('Probability $p_{1}$')
        plt.ylabel('Probability $p_{3}$')
        plt.title('Heatmap of Infectivity Variance for $p_{2} = 0.5$')
        plt.colorbar(ticks=np.arange(0, 1.6, 0.2),fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.savefig(f'Plots/Planes/Plane_Var[{pInterval}].png')

        # Determine Locations of Waves
        HighVar = VarArr[VarArr>1.0]
        with open(f'Data/WaveLocations_{pInterval}.txt', 'a') as file:
            for point in HighVar:
                Loc = np.where(VarArr==point)
                file.write(f'p1 = {pRange[Loc[1]]}, p2=[0.5], p3 = {pRange[Loc[0]]}\n')

    elif Mode == 3: # P1 CUT

        # Get Correct Interval
        Intervals = [0.01, 0.005]
        pInterval = None
        while pInterval not in Intervals:
            print(f'Available Intervals: {Intervals[0]}, {Intervals[1]}')
            pInterval = float(input('Chosen Interval: '))

        # Read in data
        Data = Table.read(f'Data/Cut_p1[{pInterval}].dat', format='ascii')
        prob = Data['Probability']
        Inf = Data['Var']
        Err = Data['Var_Err']

        # Plot cut with errorbars
        plt.figure(figsize=(8,8))
        plt.errorbar(prob, Inf, yerr=Err, fmt='none', color='black', capsize=3)
        plt.plot(prob, Inf, 'x-', lw=.5)
        plt.xlabel('Probability $p_{1}$')
        plt.ylabel('Variance $\langle I^{2} \\rangle - \langle I \\rangle ^{2}$ / $N$')
        plt.title('Plot of Infectivity Variance for $p_{2} = 0.5$, $p_{3} = 0.5$')
        plt.savefig(f'Plots/Cuts/Cut_p1[{pInterval}].png')
        plt.show()

    elif Mode == 4: # IMMUNITY FRACTION

        # Get Correct Interval
        Ns = [50, 100, 150]
        N = None
        while N not in Ns:
            print(f'Available Grid Lengths: {Ns[0]}, {Ns[1]}, {Ns[2]} (Overlay)')
            N = int(input('Chosen N: '))

        # Read in Data
        if N != 150:
            N2 = N**2
            ImmunityData = Table.read(f'Data/ImmunityData_N{N2}.dat', format='ascii')
            FracImm = ImmunityData['fImm']
            FracInf = ImmunityData['fInf']
            Error = ImmunityData['Err']
            # Determine the Zero Point
            ZeroIndex = np.where(FracInf==0)[0][0]

            # Plot
            plt.figure(figsize=(12,6))
            plt.errorbar(FracImm, FracInf, yerr=Error, fmt='none', c='black', capsize=2)
            plt.plot(FracImm, FracInf, 'x-', color='black', lw=.5, label='Fractional Infectivity')
            plt.plot(FracImm[ZeroIndex], FracInf[ZeroIndex], 'x', color='red', label=f'Minimum Immune Fraction: {FracImm[ZeroIndex].round(2)}')
            plt.legend(loc='upper right')
            plt.title(f'Plot of Fractional Infectivity as a function of Immune Fraction (for Grid Size: {N})')
            plt.xlabel('Immune Fraction $f_{Imm}$')
            plt.xlim(FracImm[0], FracImm[-1])
            plt.ylabel('Fractional Infectivity $\langle I \\rangle$ / $N$')
            plt.tight_layout()
            plt.savefig(f'Plots/Immunity/ImmuneFracPlot_N{N2}.png')

        else:
            plt.figure(figsize=(12,6))

            for N in [50, 100]:
                N2 = N**2
                if N == 50:
                    c1 = 'red'
                    c2 = 'black'
                else:
                    c1 = 'blue'
                    c2 = 'green'
                ImmunityData = Table.read(f'Data/ImmunityData_N{N2}.dat', format='ascii')
                FracImm = ImmunityData['fImm']
                FracInf = ImmunityData['fInf']
                Error = ImmunityData['Err']
                # Determine the Zero Point
                ZeroIndex = np.where(FracInf==0)[0][0]
                # Plot
                plt.errorbar(FracImm, FracInf, yerr=Error, fmt='none', c=c1, capsize=2)
                plt.plot(FracImm, FracInf, 'x-', color=c1, lw=.5, label=f'Fractional Infectivity [Grid: {N}]')
                plt.plot(FracImm[ZeroIndex], FracInf[ZeroIndex], 'x', color=c2, label=f'Minimum Immune Fraction [Grid: {N}]: {FracImm[ZeroIndex].round(2)}')

            plt.legend(loc='upper right')
            plt.title(f'Plot of Fractional Infectivity as a function of Immune Fraction (Overlay)')
            plt.xlabel('Immune Fraction $f_{Imm}$')
            plt.xlim(FracImm[0], FracImm[-1])
            plt.ylabel('Fractional Infectivity $\langle I \\rangle$ / $N$')
            plt.tight_layout()
            plt.savefig(f'Plots/Immunity/ImmuneFracPlot_Overlay.png')
