"""
This file contains all the functions used in the main program

"""
import math
import random
import numpy as np

J = 1.0

def calcdeltaE(spin,i,j):
    """
     Function to calculate delta E for a given site
     Takes in the spin configuration, and the coordinates of the site
     Returns deltaE value for the site
    """
    [lx,ly] = spin.shape
    return 2*J*spin[i,j]*(spin[(i+1)%lx, j] + spin[i,(j+1)%ly] + spin[(i-1)%lx, j] + spin[i,(j-1)%ly])


def caltotalE(spin):
    """
     Function to calculate total Energy for of a spin configuration
     Takes in spin configuration 
     Returns total energy of the spin configuration
    """
    [lx,ly] = spin.shape
    energy = 0

    for i in range(lx):
        for j in range(ly):
            Si = spin[i,j]
            Sj = spin[(i+1)%lx, j] + spin[i,(j+1)%ly] + spin[(i-1)%lx, j] + spin[i,(j-1)%ly]
            energy += -J*Si*Sj
    
    # Accounting for 4 nearest neighbours
    return energy/2


def resample(x, k, kT, N):
    """
     Function to resample a list
     Takes in the list to be resampled x, number of resamples k, kT and number of sites N
     Returns a list of c values {c1,c2,c3...ck}
    """
    c_list = []
    for i in range(k):
        y = [random.choice(x) for element in x]                    # y is a resampled list of x
        y_sq = np.square(y)                                        # y_sq is square of y
        c = (np.mean(y_sq) - np.mean(y)**2)/(N*(kT**2))            # Calculate c for the resampled list
        c_list.append(c)                                           # Add to list of c values

    # return list of c values
    return c_list

def Glauber(spin, kT, lx, ly):
    """
     Computes deltaE using the Glauber algorithm
     Takes in the spin configuration, kT, lx and ly
     Returns deltaA, and coordinates of site (itrial and jtrial)
    """

    # select spin randomly
    itrial=np.random.randint(0,lx)
    jtrial=np.random.randint(0,ly)

    # Compute delta E eg via function (account for periodic BC)            
    deltaE = calcdeltaE(spin,itrial,jtrial)

    # Return deltaE value
    return deltaE,itrial,jtrial

def Kawasaki(spin, kT, lx, ly):
    """
     Computes deltaE using the Kawasaki algorithm
     Takes in the spin configuration, kT, lx and ly
     Returns deltaA, coordinates of site A (itrialA and jtrialA) and coordinates of site B (itrialB and jtrialB)
    """
    # Selecting first random spin A
    itrialA = np.random.randint(0,50)
    jtrialA = np.random.randint(0,50)

    # Selecting second random spin B
    itrialB = np.random.randint(0,50)
    jtrialB = np.random.randint(0,50)

    # If spins are same select different B spin
    while spin[itrialA,jtrialA] == spin[itrialB,jtrialB]:
        itrialB = np.random.randint(0,50)
        jtrialB = np.random.randint(0,50)

    # Calculation of deltaE
    deltaEA = calcdeltaE(spin, itrialA, jtrialA)
    deltaEB = calcdeltaE(spin, itrialB, jtrialB)
    # Correction
    corr = 4*J*spin[itrialA,jtrialA]*spin[itrialB,jtrialB]

    # Nearest neighbour condition
    if (abs( itrialA%lx - itrialB%lx ) + abs( jtrialA%ly - jtrialB%ly )) < 2:
        deltaE = deltaEA + deltaEB - corr              
    
    # Else simple addition of Energy A and Energy B
    else:
        deltaE = deltaEA + deltaEB

    return deltaE,itrialA,jtrialA,itrialB,jtrialB