# Import Relevant Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys, time
from Lattice import Lattice
from astropy.table import Table
from matplotlib.ticker import PercentFormatter
from matplotlib import rc
rc('text', usetex=True)

##### DESCRIPTION ##############################################################

""" This method plots the histogram of stability based on the data within
    Data/HistogramData.txt                                                   """

##### MAIN METHOD ##############################################################

if __name__=="__main__":

    # Read in Histogram Data
    with open('Data/HistogramData.txt', 'r')as File:
        Data = np.array([int(line) for line in File])

    # Output Number of Data Points + Determine Maximum
    NumberOfBins = 25
    print(f'Number of DataPoints: {len(Data)}')
    HistData, Edges = np.histogram(Data, NumberOfBins)
    MinEdge = Edges[np.argmax(HistData)]
    MaxEdge = Edges[np.argmax(HistData)+1]
    Max = 0.5*(MinEdge+MaxEdge)

    # Generate Histogram Data
    plt.figure(figsize=(12,6))
    plt.axvline(Max, color='red', label=f'Est. Maximum: {int(Max)}', alpha=0.5)
    plt.hist(Data, bins=NumberOfBins, weights=np.ones(len(Data)) / len(Data), label='Binned Data')
    plt.ylabel('Probability')
    plt.xlabel('Number of Steps')
    plt.legend(loc='upper right')
    plt.title("Histogram of stabilisation steps for Conway's Game of Life")
    plt.tight_layout()
    plt.xticks(np.arange(0, 4000, 250))
    plt.xlim(0, 4000)
    plt.savefig('Plots/StabilisationHistogram.png')
