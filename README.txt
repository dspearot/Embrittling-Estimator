This code suite is a supplement for the journal article titled:
"Spectrum of Embrittling Potencies and Relation to Properties of Symmetric-Tilt Grain Boundaries"

+ This code suite assumes the density of states already exists (See "./Results/Fi_GB.csv" for the formatting).
+ The structure of the "./GBs/" folder should be preserved if additional GBs are added to the folder.
+ In the "Results/" folder of each sample there are two files. These are
	-> GBEnergies.dat 
	-> FSEnergies.dat
+ For "GBEnergies.dat" and "FSEnergies.dat" files, only first and last columns are used. Number of columns should not be modified for the script to work.
+ "GBEnergies.dat" and "FSEnergies.dat" contain the segregation energies of same atoms at a GB and a FS, respectively.
+ The order of execution:
	1) Population.py
	2) Samples.py
+ This code suite requires numpy and pandas libraries.	
