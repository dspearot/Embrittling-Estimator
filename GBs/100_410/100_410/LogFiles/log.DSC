shell cd ..

include ${GBName}.config
include 100_410.config
# Interatomic Potential
if "${IntPotVar} == 1" then  "variable 		IntPot index 'NiS.eam.alloy'"  else  "variable 		IntPot index 'NiS.ff'"
variable 		IntPot index 'NiS.eam.alloy'

variable 		a0 equal 3.52
variable 		ebulkatom equal -4.44999999950988

# Lattice orientations
variable		gb1x1 equal 1
variable		gb1x2 equal -4
variable		gb1x3 equal 0
variable		gb1y1 equal 4
variable		gb1y2 equal 1
variable		gb1y3 equal 0
variable		gb1z1 equal 0
variable		gb1z2 equal 0
variable		gb1z3 equal 1
variable		gb2x1 equal -1
variable		gb2x2 equal 4
variable		gb2x3 equal 0
variable		gb2y1 equal 4
variable		gb2y2 equal 1
variable		gb2y3 equal 0
variable		gb2z1 equal 0
variable		gb2z2 equal 0
variable		gb2z3 equal -1

# Variables for S insertion
variable        deltamu equal 0.1 #eV
variable        sigma string 410
variable        GBRegion1 equal 53.1249
variable        GBHalfRegionUpper equal 110.23
variable 		GBHalfRegionLower equal -2.37
variable        slice equal 5.0
variable        t equal 948
variable        seed equal 63241
#variable 		con equal 7.86 #S/nm^2
#variable        mu equal -0.4 #eV
#variable        bulkS equal -2.9371861774248611

# Find the repeating distances
variable 		square1 equal ${gb1x1}^2+${gb1x2}^2+${gb1x3}^2
variable 		square1 equal 1^2+${gb1x2}^2+${gb1x3}^2
variable 		square1 equal 1^2+-4^2+${gb1x3}^2
variable 		square1 equal 1^2+-4^2+0^2
variable 		square2 equal ${gb1y1}^2+${gb1y2}^2+${gb1y3}^2
variable 		square2 equal 4^2+${gb1y2}^2+${gb1y3}^2
variable 		square2 equal 4^2+1^2+${gb1y3}^2
variable 		square2 equal 4^2+1^2+0^2
variable 		square3 equal ${gb1z1}^2+${gb1z2}^2+${gb1z3}^2
variable 		square3 equal 0^2+${gb1z2}^2+${gb1z3}^2
variable 		square3 equal 0^2+0^2+${gb1z3}^2
variable 		square3 equal 0^2+0^2+1^2

variable repDistx equal 4.123105625617661
variable repDisty equal 4.123105625617661
variable repDistz equal 1.0

# 3 Layer Structure dimensions (Step1 & Step2)
variable		xlo3Layer equal -0.001
variable		xhi3Layer equal 1*${a0}*${repDistx}
variable		xhi3Layer equal 1*3.52*${repDistx}
variable		xhi3Layer equal 1*3.52*4.12310562561766
variable		ylo3Layer equal -0.001
variable		yhi3Layer equal 6*${a0}*${repDisty}
variable		yhi3Layer equal 6*3.52*${repDisty}
variable		yhi3Layer equal 6*3.52*4.12310562561766
variable		zlo3Layer equal -0.001
variable		zhi3Layer equal 1*${a0}*${repDistz}
variable		zhi3Layer equal 1*3.52*${repDistz}
variable		zhi3Layer equal 1*3.52*1
variable        GB1 equal ${yhi3Layer}/4
variable        GB1 equal 87.079990813045/4
variable        GB2 equal ${yhi3Layer}/4*3
variable        GB2 equal 87.079990813045/4*3

# Iteration Numbers
variable		iter_x equal 10
variable		iter_y equal 10
variable		iter_z equal 10

variable		raster_x equal 10

# Box Dimensions for in.Structures and in.SGCMC
variable		xlo equal -0.001
variable		xhi equal 1*${a0}*${repDistx}
variable		xhi equal 1*3.52*${repDistx}
variable		xhi equal 1*3.52*4.12310562561766
variable		ylo equal -0.001
variable		yhi equal 4*${a0}*${repDisty}
variable		yhi equal 4*3.52*${repDisty}
variable		yhi equal 4*3.52*4.12310562561766
variable		zlo equal -0.001
variable		zhi equal 1*${a0}*${repDistz}
variable		zhi equal 1*3.52*${repDistz}
variable		zhi equal 1*3.52*1
variable		GB1Str equal 2*${a0}*${repDisty}
variable		GB1Str equal 2*3.52*${repDisty}
variable		GB1Str equal 2*3.52*4.12310562561766
variable		FS1Str equal 3*${a0}*${repDisty}
variable		FS1Str equal 3*3.52*${repDisty}
variable		FS1Str equal 3*3.52*4.12310562561766

shell cd DumpFiles
shell rm dump.DSC.*
WARNING: Shell command 'rm' failed with error 'No such file or directory' (../input.cpp:1313)
shell cd ..

# ---------- Initialize Simulation ---------------------
clear

if "${IntPotVar} == 1" then "units metal" "atom_style atomic" else "units real" "atom_style charge"
units metal
atom_style atomic

dimension 3
boundary p p p

# ---------- Variables ----------------------------------
variable		xlo equal -0.001
variable		xhi equal 2*${a0}*${repDistx}
variable		xhi equal 2*3.52*${repDistx}
variable		xhi equal 2*3.52*4.12310562561766
variable		ylo equal -0.001
variable		yhi equal 2*${a0}*${repDisty}
variable		yhi equal 2*3.52*${repDisty}
variable		yhi equal 2*3.52*4.12310562561766
variable		zlo equal -${a0}*${repDistz}/4
variable		zlo equal -3.52*${repDistz}/4
variable		zlo equal -3.52*1/4
variable		zhi equal ${a0}*${repDistz}/4
variable		zhi equal 3.52*${repDistz}/4
variable		zhi equal 3.52*1/4

variable		boxlengthx equal ${xhi}
variable		boxlengthx equal 29.0266636043483
variable		boxlengthy equal ${yhi}
variable		boxlengthy equal 29.0266636043483
variable		boxlengthz equal ${zhi}
variable		boxlengthz equal 0.88

# ---------- Create Atomistic Structure ---------------------
lattice		fcc ${a0}
lattice		fcc 3.52
Lattice spacing in x,y,z = 3.52 3.52 3.52
region		whole block ${xlo} ${boxlengthx} ${ylo} ${boxlengthy} ${zlo} ${boxlengthz} units box
region		whole block -0.001 ${boxlengthx} ${ylo} ${boxlengthy} ${zlo} ${boxlengthz} units box
region		whole block -0.001 29.0266636043483 ${ylo} ${boxlengthy} ${zlo} ${boxlengthz} units box
region		whole block -0.001 29.0266636043483 -0.001 ${boxlengthy} ${zlo} ${boxlengthz} units box
region		whole block -0.001 29.0266636043483 -0.001 29.0266636043483 ${zlo} ${boxlengthz} units box
region		whole block -0.001 29.0266636043483 -0.001 29.0266636043483 -0.88 ${boxlengthz} units box
region		whole block -0.001 29.0266636043483 -0.001 29.0266636043483 -0.88 0.88 units box
create_box	2 whole
Created orthogonal box = (-0.001 -0.001 -0.88) to (29.0267 29.0267 0.88)
  4 by 8 by 1 MPI processor grid

# Lattice 1
lattice			fcc ${a0} orient x ${gb1x1} ${gb1x2} ${gb1x3} orient y ${gb1y1} ${gb1y2} ${gb1y3} orient z ${gb1z1} ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x ${gb1x1} ${gb1x2} ${gb1x3} orient y ${gb1y1} ${gb1y2} ${gb1y3} orient z ${gb1z1} ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x 1 ${gb1x2} ${gb1x3} orient y ${gb1y1} ${gb1y2} ${gb1y3} orient z ${gb1z1} ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x 1 -4 ${gb1x3} orient y ${gb1y1} ${gb1y2} ${gb1y3} orient z ${gb1z1} ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x 1 -4 0 orient y ${gb1y1} ${gb1y2} ${gb1y3} orient z ${gb1z1} ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x 1 -4 0 orient y 4 ${gb1y2} ${gb1y3} orient z ${gb1z1} ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x 1 -4 0 orient y 4 1 ${gb1y3} orient z ${gb1z1} ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x 1 -4 0 orient y 4 1 0 orient z ${gb1z1} ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x 1 -4 0 orient y 4 1 0 orient z 0 ${gb1z2} ${gb1z3}
lattice			fcc 3.52 orient x 1 -4 0 orient y 4 1 0 orient z 0 0 ${gb1z3}
lattice			fcc 3.52 orient x 1 -4 0 orient y 4 1 0 orient z 0 0 1
Lattice spacing in x,y,z = 4.26863 4.26863 3.52
region		lower block ${xlo} ${xhi} ${ylo} ${yhi} ${zlo} ${zhi} units box
region		lower block -0.001 ${xhi} ${ylo} ${yhi} ${zlo} ${zhi} units box
region		lower block -0.001 29.0266636043483 ${ylo} ${yhi} ${zlo} ${zhi} units box
region		lower block -0.001 29.0266636043483 -0.001 ${yhi} ${zlo} ${zhi} units box
region		lower block -0.001 29.0266636043483 -0.001 29.0266636043483 ${zlo} ${zhi} units box
region		lower block -0.001 29.0266636043483 -0.001 29.0266636043483 -0.88 ${zhi} units box
region		lower block -0.001 29.0266636043483 -0.001 29.0266636043483 -0.88 0.88 units box
create_atoms	1 region lower
Created 136 atoms
  create_atoms CPU = 0.00154929 secs

# Lattice 2
lattice			fcc ${a0} orient x ${gb2x1} ${gb2x2} ${gb2x3} orient y ${gb2y1} ${gb2y2} ${gb2y3} orient z ${gb2z1} ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x ${gb2x1} ${gb2x2} ${gb2x3} orient y ${gb2y1} ${gb2y2} ${gb2y3} orient z ${gb2z1} ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x -1 ${gb2x2} ${gb2x3} orient y ${gb2y1} ${gb2y2} ${gb2y3} orient z ${gb2z1} ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x -1 4 ${gb2x3} orient y ${gb2y1} ${gb2y2} ${gb2y3} orient z ${gb2z1} ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x -1 4 0 orient y ${gb2y1} ${gb2y2} ${gb2y3} orient z ${gb2z1} ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x -1 4 0 orient y 4 ${gb2y2} ${gb2y3} orient z ${gb2z1} ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x -1 4 0 orient y 4 1 ${gb2y3} orient z ${gb2z1} ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x -1 4 0 orient y 4 1 0 orient z ${gb2z1} ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x -1 4 0 orient y 4 1 0 orient z 0 ${gb2z2} ${gb2z3}
lattice			fcc 3.52 orient x -1 4 0 orient y 4 1 0 orient z 0 0 ${gb2z3}
lattice			fcc 3.52 orient x -1 4 0 orient y 4 1 0 orient z 0 0 -1
Lattice spacing in x,y,z = 4.26863 4.26863 3.52
region		middle block ${xlo} ${xhi} ${ylo} ${yhi} ${zlo} ${zhi} units box
region		middle block -0.001 ${xhi} ${ylo} ${yhi} ${zlo} ${zhi} units box
region		middle block -0.001 29.0266636043483 ${ylo} ${yhi} ${zlo} ${zhi} units box
region		middle block -0.001 29.0266636043483 -0.001 ${yhi} ${zlo} ${zhi} units box
region		middle block -0.001 29.0266636043483 -0.001 29.0266636043483 ${zlo} ${zhi} units box
region		middle block -0.001 29.0266636043483 -0.001 29.0266636043483 -0.88 ${zhi} units box
region		middle block -0.001 29.0266636043483 -0.001 29.0266636043483 -0.88 0.88 units box
create_atoms	2 region middle
Created 136 atoms
  create_atoms CPU = 7.15998e-05 secs

# ---------- Define Interatomic Potential ---------------------
if "${IntPotVar} == 1" then "pair_style      eam/alloy" else "pair_style      reax/c NULL checkqeq no"
pair_style      eam/alloy

pair_coeff * * ../../${IntPot} Ni Ni
pair_coeff * * ../../NiS.eam.alloy Ni Ni

if "${IntPotVar} == 0" then "mass 1 58.71" "mass 2 58.71"

shell cd DumpFiles
dump    1 all custom 1 dump.DSC.* id type x y z
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Neighbor list info ...
  update every 1 steps, delay 10 steps, check yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 8.5
  ghost atom cutoff = 8.5
  binsize = 4.25, bins = 7 7 1
  1 neighbor lists, perpetual/occasional/extra = 1 0 0
  (1) pair eam/alloy, perpetual
      attributes: half, newton on
      pair build: half/bin/atomonly/newton
      stencil: half/bin/3d/newton
      bin: standard
WARNING: Proc sub-domain size < neighbor skin, could lead to lost atoms (../domain.cpp:931)
Per MPI rank memory allocation (min/avg/max) = 4.506 | 4.506 | 4.506 Mbytes
Step Temp E_pair E_mol TotEng Press 
       0            0         -nan            0         -nan         -nan 
Loop time of 2.11612e-05 on 32 procs for 0 steps with 272 atoms

100.4% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 2.116e-05  |            |       |100.00

Nlocal:    8.5 ave 10 max 8 min
Histogram: 24 0 0 0 0 0 0 0 0 8
Nghost:    1729.5 ave 1730 max 1728 min
Histogram: 8 0 0 0 0 0 0 0 0 24
Neighs:    1964.25 ave 2317 max 1840 min
Histogram: 24 0 0 0 0 0 0 0 4 4

Total # of neighbors = 62856
Ave neighs/atom = 231.088
Neighbor list builds = 0
Dangerous builds = 0
undump 1
shell cd ..
Total wall time: 0:00:00
