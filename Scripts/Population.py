# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# This code is a supplement for the journal article titled:
#   "Spectrum of Embrittling Potencies and Relation to Properties of
#   Symmetric-Tilt Grain Boundaries"
# ------------------
# This code performs the following tasks:
# 1) Obtains density of states from the previous step
# 2) Calculates Xi and Pi (check the paper for definitions) at the population
#  level (Fig.4)
# 3) Write Xi and Pi calculated in this step to a data frame, to be processed
#  at the sample level
# --- Definitions and Abbreviations --
# GB: Grain boundary
# FS: Free surface
# ------------------
# Authors: Doruk Aksoy (1), Rémi Dingreville (2), Douglas E. Spearot (1,*)
# (1) University of Florida, Gainesville, FL, USA
# (2) Center for Integrated Nanotechnologies, Sandia National Laboratories,
#  Albuquerque, NM, USA
# (*) dspearot@ufl.edu
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#%% Imports
import numpy as np
import pandas as pd

# %% Define functions
def calcXtot(delta_E_seg_GB_i,Fi,X_bulk):
    '''
    Calculate total solute concentration from bulk solute concentration.

    Parameters
    ----------
    X_bulk : Bulk solute concentration
    delta_E_seg_GB_i : All segregation energies for each site type i
    Fi : Density of states for each site type within the population

    Returns
    -------
    X_tot : Total solute concentration within the population

    '''
    # Number of site types
    n_site_types = np.size(Fi,axis=0)
    # Initialize and calculate the probability distribution function for each
    # site type i with the given bulk concentration
    Xi_with_bulk = np.zeros(n_site_types)
    for i in range(n_site_types): Xi_with_bulk[i] = 1 / (1 + ((1 - X_bulk) / X_bulk) * np.exp( delta_E_seg_GB_i[i] / (kB * T)))
    # Calculate the effective solute concentration
    X_bar = np.sum(Fi * Xi_with_bulk)
    # Return the total solute concentration
    return ((1 - f_int) * X_bulk + f_int * X_bar)

def fromXtotToXbulk(delta_E_seg_GB_i,Fi,X_tot,tol):
    '''
    Calculate bulk solute concentration from total solute concentration using
    midpoint trial and improvement solver.

    Parameters
    ----------
    delta_E_seg_GB_i : All segregation energies for each site type i
    Fi : Density of states for each site type
    X_tot : Total solute concentration
    tol : Tolerance

    Returns
    -------
    If a result is found, return X_bulk.

    '''
    # Initial lower and upper estimates
    x_lo = 0.0
    x_hi = X_tot*2
    # Initial guess
    x_0 = (x_lo + x_hi)/2
    # Calculate a trial value using calcXtot function
    X_tot_trial = calcXtot(delta_E_seg_GB_i,Fi,x_0)
    # Initialize iteration counter
    iter_count = 0
    # Maximum number of iterations
    max_iter = 100
    # Check if the result is within the tolerance and number of iterations
    # is less than the maximum value
    while((np.abs(X_tot_trial - X_tot) > tol) and (iter_count < max_iter)):
        if(X_tot_trial > X_tot):
            x_hi = x_0
            x_0 = (x_hi + x_lo)/2 # Next guess
        else:
            x_lo = x_0
            x_0 = (x_hi + x_lo)/2 # Next guess
        # Calculate the new trial value using calcXtot function
        X_tot_trial = calcXtot(delta_E_seg_GB_i,Fi,x_0)
        # Increment the iteration counter
        iter_count +=1
    # Check whether a total solute concentration can be found
    if (iter_count == max_iter):
        print("Could not find a value.")
        return (0)
    else:
        return (x_0)

def calcPopProp(delta_E_seg_GB_i,Fi,X_tot):
    '''
    Calculate population properties.

    Parameters
    ----------
    delta_E_seg_GB_i : All segregation energies for each site type i
    Fi : Density of states for each site type
    X_tot : Total solute concentration

    Returns
    -------
    X_bulk : Bulk solute concentration
    Xi : Fraction of occupied type i sites
    Pi : Solute occupancy density
    X_bar : Effective solute concentration
    delta_E_bar_seg_GB_i : Effective segregation energy per site type i
    delta_E_bar_seg_GB : Total effective segregation energy

    '''
    # Calculate the bulk concentration using fromXtotToXbulk function
    X_bulk = fromXtotToXbulk(delta_E_seg_GB_i,Fi,X_tot,1E-4)
    # Raise an exception if a bulk solute concentration cannot be calculated with given total solute concentration
    if (X_bulk==0):
        raise Exception('Error: Cannot calculate a bulk solute concentration with given total solute concentration.')
    # Calculate the site specific probability distribution function and convert it to numpy array
    Xi = [(1/(1+ ((1-X_bulk)/X_bulk) * np.exp( delta_E_seg_GB_i[i] / (kB*T)))) for i in range(np.size(delta_E_seg_GB_i))]
    Xi = np.array(Xi)
    # Site occupancy
    Pi = Fi * Xi
    # Effective solute concentration
    X_bar = np.sum(Pi)
    # Effective segregation energy for each site type i
    delta_E_bar_seg_GB_i = (1/(X_bar*(1-X_bar))) * (Fi * delta_E_seg_GB_i * Xi * (1-Xi))
    # Effective segregation energy
    delta_E_bar_seg_GB =  np.sum(delta_E_bar_seg_GB_i)
    # Return all calculated properties
    return (X_bulk,Xi,Pi,X_bar,delta_E_bar_seg_GB_i,delta_E_bar_seg_GB)

# %% MAIN
# Read-in normalized density of states (Format: Index/Energies/Frequencies)
df_Fi_GB = pd.read_csv("../Results/Fi_GB.csv",index_col = 0)

# Segregation energies for each site type i
delta_E_seg_GB_i = np.array(df_Fi_GB['Energy'])
# Density of states
Fi = np.array(df_Fi_GB['Freq'])

# %% Variables
# Total solute concentration
X_tot = 15/100 # no of solute atoms/no of GB atoms
# Fraction of interface sites to all segregation sites
f_int = 0.162
# Boltzmann Constant in eV K-1
kB = 0.00008617333262
# Temperature
T = 300 # K

# %% Calculate properties corresponding to the GB population using calcPopProp function
(X_bulk,Xi,Pi,X_bar,delta_E_bar_seg_GB_i,delta_E_bar_seg_GB) = calcPopProp(delta_E_seg_GB_i,Fi,X_tot)

# %% Create a data frame with the population properties
df_Pop = pd.DataFrame(np.transpose([delta_E_seg_GB_i, Fi, Xi, Pi]),columns=['delta_E_seg_GB_i','Fi','Xi','Pi']).astype(float)
# Convert data frame to csv
df_Pop.to_csv("../Results/Pop.csv")
