"""
Program to plot the various graphs of observables vs temperature
Uses the Observables.txt file to read columns of data and plots them
Column order : Temperature, Average Energy, Average |M|, Specific heat, Error in Specific Heat, Susceptibility

"""

import numpy as np
import matplotlib.pyplot as pyplot

# Load file
data = np.loadtxt("Observables.txt")

# Plot observables vs Temperature

# Mean Energy
pyplot.title('Average Energy vs. Temperature')
pyplot.xlabel('Temperature')
pyplot.ylabel('Energy')
pyplot.plot(data[:,0], data[:,1], linestyle='--', markeredgecolor = 'k', marker='x', color='b')
pyplot.show()

# Mean |M|
pyplot.title('Average |M| vs. Temperature')
pyplot.xlabel('Temperature')
pyplot.ylabel('|M|')
pyplot.plot(data[:,0], data[:,2], linestyle='--', markeredgecolor = 'k', marker='x', color='b')
pyplot.show()

# Specific Heat with error bars
fig = pyplot.figure()
pyplot.title('C vs. Temperature')
pyplot.xlabel('Temperature')
pyplot.ylabel('Specific Heat')
pyplot.errorbar(data[:,0], data[:,3],yerr = data[:,4], markeredgecolor = 'k', ecolor = 'r', capsize = 2, linestyle='--', elinewidth = 0.5, marker='x', color='b')
pyplot.show()

# Susceptibility
pyplot.title('Susceptibility vs. Temperature')
pyplot.xlabel('Temperature')
pyplot.ylabel('Susceptibility')
pyplot.plot(data[:,0], data[:,5], linestyle='--', marker='x', markeredgecolor = 'k', color='b')
pyplot.show()