"""
Program to simulate the Ising model for a range of temperature values using either Glauber or Kawasaki dynamics
Takes in length N of lattice where the lattice is a square lattice of size N X N
Takes in value for start temperature, end temperature and number of equidistant temperature values in between
kB and J are in reduced units = 1.0

Prompts from user 
1) For Glauber dynamics
2) For kawasaki

Performs a sweep of N X N calculations and plots the visualisation every 10 sweeps
Equilibrium time is set to 200 sweeps
Observables total energy of configration, magnetisation, absolute magnetisation are calculated every 10 sweeps

Susceptibility and specific heat are calculated after sweeps are completed
Error in specific heat is calculated using bootstrap resampling method

Temperature, average energy, average |M|, specific heat, error in specific heat and susceptibility are written to file Observables.txt for each temperature iteration


"""


from typing import no_type_check
import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from methods import calcdeltaE, caltotalE, resample, Glauber, Kawasaki


##############   MAIN FUNCTION   #######################

J = 1.0
nstep = 1000              # Number of sweeps

# Input

if(len(sys.argv) != 5):
   print ("Usage python CP1_loop.py N T_start T_end number_of_temps")
   sys.exit()

lx=int(sys.argv[1]) 
ly=lx 

# Number of sites
N = lx*ly

# Dynamics = 1 for Glauber, 2 for Kawasaki
dynamics = int(input("Select dynamics\n1)Glauber\n2)Kawasaki\n"))
if dynamics != 1 and dynamics != 2:
    print("Must be 1 or 2")
    sys.exit()

# Set up files
mfile = open("Observables.txt","w")

# To calculate time elapsed
start = time.time()

# Creating list of temperature values
temp_low = float(sys.argv[2])                                   # start temp
temp_high = float(sys.argv[3])                                  # end temp
no_of_temps = int(sys.argv[4])                                  # Number of equidistant temp values 
Temp_list = np.linspace(temp_low,temp_high,num = no_of_temps)   # List of temperatures to iterate over
Temp_list = np.round(Temp_list,3)                               # limiting to three decimal places

# Loop through all temperature values
for kT in Temp_list:

    # Initialise lists to store observables
    mag_list = []                               # List for magnetisation values
    magsq_list = []                             # List for magnetisation squared values
    energy_list = []                            # List for energy values
    energysq_list = []                          # List for energy squared values

    # Initialise spins randomly
    spin=np.zeros((lx,ly),dtype=float)
    
    for i in range(lx):
        for j in range(ly):
            r=random.random()
            if(r<0.5):
                spin[i,j] = -1
            else:
                spin[i,j] = 1

    # Set up plot
    fig = plt.figure()
    im = plt.imshow(spin, animated=True)

    # Start sweeps            
    for n in range(nstep):

        # Looping through lx*ly iterations
        for i in range(lx):
            for j in range(ly):

                # Glauber
                if (dynamics == 1):
                    # Get deltaE and coordinates of site from Glauber function
                    deltaE,itrial,jtrial = Glauber(spin, kT, lx, ly)   

                    # Perform metropolis test
                    r = random.random()
                    # If random number is less than exp(-deltaE/kT), switch spin
                    if(r<+np.exp(-deltaE/kT)):
                        spin[itrial,jtrial] = -spin[itrial,jtrial]

                # Kawasaki
                else:
                    # Get deltaE and coordinates of sites A and B from function
                    deltaE,itrialA,jtrialA,itrialB,jtrialB = Kawasaki(spin, kT, lx, ly) 

                    # Perform metropolis test
                    r = random.random()
                    # If random number is less than exp(-deltaE/kT), switch spin
                    if(r<+np.exp(-deltaE/kT)):
                        spin[itrialA,jtrialA] = -spin[itrialA,jtrialA]
                        spin[itrialB,jtrialB] = -spin[itrialB,jtrialB]
                    
        # Occasionally plot or update measurements, eg every 10 sweeps (Auto correlation time)
        if(n%10==0): 
            print(kT,n)
            # Plot
            plt.cla()
            im=plt.imshow(spin, animated=True)
            plt.draw()
            plt.pause(0.0001)

            # Wait 200 sweeps for equilibrium
            if n > 200:
                # Magnetisation
                M = np.sum(spin)
                mag_list.append(M)
                magsq_list.append(M**2)
                # Total Energy of spin configuration
                totalE = caltotalE(spin)
                energy_list.append(totalE)
                energysq_list.append(totalE**2)
                


    mag_abs_list = [abs(element) for element in mag_list]     # List of absolute values of magnetisation
    exp_Mabs = np.mean(mag_abs_list)                          # Expected value of absolute value of magnetisation
    exp_M = np.mean(mag_list)                                 # Expected value of magnetisation
    exp_Msq = np.mean(magsq_list)                             # Expected value of magnetisation square 
    exp_E = np.mean(energy_list)                              # Expected value of energy
    exp_Esq = np.mean(energysq_list)                          # Expected value of energy square

    # Susceptibility
    susc = (exp_Msq-exp_M**2)/(N*kT)
    # Specific Heat
    C = (exp_Esq-exp_E**2)/(N*(kT**2))

    ###### Error in specific heat (Resampling using bootstrap method) ################
    C_list = resample(energy_list, 100, kT, N)
    err = np.std(C_list)

    # Write observables to file
    mfile.write('%.1f %.1f %d %f %f %f \n' %(kT,exp_E,exp_Mabs,C,err,susc))


print("Time Elapsed = {:.2f}s".format(time.time()-start))

# Close files
mfile.close()
