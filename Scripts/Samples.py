# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# This code is a supplement for the journal article titled:
#   "Spectrum of Embrittling Potencies and Relation to Properties of
#   Symmetric-Tilt Grain Boundaries"
# ------------------
# This code performs the following tasks:
# 1) Reads in Fi, Xi, Pi from the previous step
# 2) Calculates site-specific properties that are shown in Table 2 and Fig. 6
# 3) Calculates collective-behavior properties that are shown in Table 3 and Fig. 5
# 4) Generates all data frames for plotting
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
from os import listdir,path

# %% Define functions
def getNumOfAtoms(file_path, file_name):
    '''
    Obtain number of atoms from the file.

    Parameters
    ----------
    file_path : File path
    file_name : Name of the file

    Returns
    -------
    Number of atoms

    '''
    with open(path.join(file_path,file_name), 'r') as atoms_file:
        # Number of atoms is equal to number of lines without the header
        lineCount = 0
        for line in atoms_file:
            lineCount += 1
    return int(lineCount)-1

def getEnergies(file_path, file_name, arr):
    '''
    Function to obtain energies from file

    Parameters
    ----------
    file_path : File path
    file_name : Name of the file
    arr : Array to write energies

    '''
    with open(path.join(file_path,file_name), 'r') as results_file:
        for ind,line in enumerate(results_file):
            # Skip the header
            if "#" not in line:
                line = line.split()
                for j in range(int(np.size(line))):
                    arr[int(ind)-1,j] = line[j]

def segEngOcc(energies,col_num,NDIGITS):
    '''
    Function to obtain energies from file

    Parameters
    ----------
    energies : Energies obtained from simulations
    col_num : Segregation energy column number
    NDIGITS : Number of digits to consider when looking at unique segregation
        energies

    Returns
    -------
    DE_seg_i_GB : Segregation energy of site type i
    N_hat_i_GB : Number of occurences of the segregation energy of site type i
    num_site_types : Total number of unique site types
    site_type_ind : Indices of matching energies between DE_seg_i_GB array and energies array

    '''

    # Round energies by the given number of digits, and then find number of unique energies and number of occurences
    DE_seg_i_GB,N_hat_i_GB = np.unique(np.round(energies[np.nonzero(energies[:,col_num]),col_num],NDIGITS), return_counts=True)
    # Number of site types
    num_site_types = int(np.size(DE_seg_i_GB))
    # We will use the site_type_ind list to match the site types between GBs and FSs.
    site_type_ind = []
    # Now that we have matched the rounded energies, find originals and put back into DE_seg_i_GB array
    for i in range(num_site_types):
        site_type_ind.append(np.where(np.round(energies[np.nonzero(energies[:,col_num]),col_num],NDIGITS) == DE_seg_i_GB[i])[1][0])
        DE_seg_i_GB[i] = energies[site_type_ind[i],col_num]
    return (DE_seg_i_GB, N_hat_i_GB, num_site_types, site_type_ind)
# %% MAIN
# Read in data frames
df_Pop = pd.read_csv("../Results/Pop.csv",index_col = 0).astype(float)

# From data frame to arrays
delta_E_seg_GB_i = np.array(df_Pop['delta_E_seg_GB_i'])
Pi = np.array(df_Pop['Pi'])

# Round by this number when comparing energies
NDIGITS = 3

# Perform simulations for all given models
allSims = listdir('../GBs/')

# %% Create a data frame to store all results
# Define columns (three properties shown in Fig. 5)
columns_all = ["DE_hat_b","PR_hat_GB","E_hat_b"]
# Tilt and GB normals as indices of the data frame
tilt_axes = [sim.split('_')[0] for sim in allSims]
GB_normals = [sim.split('_')[1] for sim in allSims]
# Levels required for a multi index data frame
levels_all = list(zip(*[tilt_axes, GB_normals]))
# Define indices
index_all = pd.MultiIndex.from_tuples(levels_all, names=['Tilt', 'Normal'])
# Initialize the data frame
df_all = pd.DataFrame(index = index_all, columns=columns_all)

#%% For each sample
for indSim,sim in enumerate(allSims):

    # Obtain GB normal and tilt axes from the folder names
    GB_normal = str(sim.split('_')[1])
    GB_tilt = str(sim.split('_')[0])

    # Model path
    model_path = path.join("../GBs/", str(sim) + "/")

    # Read in number of GB atoms considered in the simulation
    N_hat_GB = getNumOfAtoms(path.join(model_path, "Results/"),"GBEnergies.dat")

    # Initialize an array for energies of individual sites in GB models
    GBenergies = np.zeros((N_hat_GB,5))
    # Initialize an array for energies of individual sites in FS models
    FSenergies = np.zeros((N_hat_GB,5))

    try:
        # Read energies for each sample
        getEnergies(path.join(model_path, "Results/"),"GBEnergies.dat",GBenergies)
        getEnergies(path.join(model_path, "Results/"),"FSEnergies.dat",FSenergies)

        # Sort by atom ID
        GBenergies = GBenergies[np.argsort(GBenergies[:,0]),:]
        FSenergies = FSenergies[np.argsort(FSenergies[:,0]),:]

        # Weed out non-matching simulations (if one of two simulations per atom ID is failed)
        # Find out the intersection vector of two arrays, then delete rows with different atom IDs
        for ind,val in enumerate(np.asarray(np.intersect1d(GBenergies[:,0],FSenergies[:,0]),dtype=int)):
            if (not np.searchsorted(GBenergies[:,0],val) == ind):
                GBenergies = np.delete(GBenergies,ind,0)
            if (not np.searchsorted(FSenergies[:,0],val) == ind):
                FSenergies = np.delete(FSenergies,ind,0)

        # Update number of atoms
        N_hat_GB = np.size(GBenergies,axis=0)

        # Find unique segregation energies and their number of occurences using segEngOcc function
        DE_seg_i_GB, N_hat_i_GB, num_site_types_GB, site_type_ind_GB  = segEngOcc(GBenergies,4,NDIGITS)
        # Site type indices should be preserved after cleavage (See Section 4)
        DE_seg_i_FS = FSenergies[site_type_ind_GB,4]
        # Embrittling potencies
        DE_b_i = GBenergies[site_type_ind_GB,4]-FSenergies[site_type_ind_GB,4]

        # Site occupancies
        P_bar_i_GB = np.zeros(num_site_types_GB)

        # Obtain P_bar_i_GB from the population (closest value)
        for i in range(num_site_types_GB): P_bar_i_GB[i] = Pi[(np.abs(delta_E_seg_GB_i - DE_seg_i_GB[i])).argmin()]

        # Rescaled site occupancy for each site type i
        PR_hat_i_GB = P_bar_i_GB/np.sum(np.multiply(P_bar_i_GB, N_hat_i_GB))

        # Site specific embrittling estimator
        E_hat_b_i = np.multiply(PR_hat_i_GB,DE_b_i)

        # Sample embrittling estimator
        E_hat_b = np.sum(np.multiply(np.multiply(PR_hat_i_GB,N_hat_i_GB),DE_b_i))/(N_hat_GB)

        # Write properties to the all results data frame
        df_all['DE_hat_b'][GB_tilt,GB_normal] = np.sum(np.mean(np.multiply(DE_b_i,N_hat_i_GB)))/N_hat_GB
        df_all['PR_hat_GB'][GB_tilt,GB_normal] = np.sum(np.mean(np.multiply(PR_hat_i_GB,N_hat_i_GB)))/N_hat_GB
        df_all['E_hat_b'][GB_tilt,GB_normal] = E_hat_b

    except:
        print(indSim+1,sim,"Properties not calculated!")
        continue

# %% To csv
df_all.to_csv("../Results/AllResults.csv")