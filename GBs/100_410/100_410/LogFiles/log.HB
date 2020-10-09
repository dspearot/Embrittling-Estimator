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

#####################################################
#################       (f)       ###################
#####################################################
#log 	log.GBBulk

# --------------- INITIALIZATION ------------------
clear
if "${IntPotVar} == 1" then "units metal" "atom_style atomic" else "units real" "atom_style charge"
units metal
atom_style atomic

dimension 3
boundary p p p

shell cd RestartFiles
read_data 	data.GBBulk
  orthogonal box = (0.004876 0.034254 0.000425) to (29.01 87.0447 7.03672)
  4 by 8 by 1 MPI processor grid
  reading atoms ...
  1616 atoms
  read_data CPU = 0.00345102 secs
shell cd ..
# ---------- Define Interatomic Potential ---------------------
if "${IntPotVar} == 1" then "pair_style      eam/alloy" else "pair_style      reax/c NULL safezone 1.6 mincap 100"
pair_style      eam/alloy

pair_coeff * * ../../${IntPot} Ni S Ni Ni
pair_coeff * * ../../NiS.eam.alloy Ni S Ni Ni

if "${IntPotVar} == 0" then "mass 1 58.71" "mass 2 32.06" "mass 1 58.71" "mass 1 58.71"

# ---------- Groups ---------------------------------
group 			nickel type 1
1608 atoms in group nickel
group 			sulfur type 2
8 atoms in group sulfur

#--------------------------- Computations ----------------------------
compute csym all centro/atom fcc
compute PE all pe/atom
compute eatoms all reduce sum c_PE

if "${IntPotVar} == 0" then "fix QEQREAX all qeq/reax 1 0.0 10.0 1.0e-6 ../../paramNiS.qeq"

# ------------ Thermo style ---------------------------------
thermo 			100
thermo_style    custom step temp pe c_eatoms
thermo_modify   format float %24.16g

# ----------- Dump Results ---------------------------------
shell cd DumpFiles
dump            GBBulkDump all custom 100 dump.GBBulk id type x y z c_PE c_csym
variable NVTVar equal 0
variable BRVar equal 1
include ../../../min.config
min_style cg
minimize 1e-25 1e-25 1000000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Neighbor list info ...
  update every 1 steps, delay 0 steps, check yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 8.5
  ghost atom cutoff = 8.5
  binsize = 4.25, bins = 7 21 2
  2 neighbor lists, perpetual/occasional/extra = 1 1 0
  (1) pair eam/alloy, perpetual
      attributes: half, newton on
      pair build: half/bin/atomonly/newton
      stencil: half/bin/3d/newton
      bin: standard
  (2) compute centro/atom, occasional
      attributes: full, newton on
      pair build: full/bin/atomonly
      stencil: full/bin/3d
      bin: standard
Per MPI rank memory allocation (min/avg/max) = 6.142 | 6.143 | 6.144 Mbytes
Step Temp PotEng c_eatoms 
       0                        0        -7146.28377393132       -7146.283773931321 
     100                        0       -7146.523648686853       -7146.523648686852 
     135                        0       -7146.523648687695       -7146.523648687694 
Loop time of 0.189953 on 32 procs for 135 steps with 1616 atoms

99.9% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -7146.28377393     -7146.52364869     -7146.52364869
  Force two-norm initial, final = 2.34985 3.97794e-06
  Force max component initial, final = 0.306536 5.36243e-07
  Final line search alpha, max atom move = 1 5.36243e-07
  Iterations, force evaluations = 135 270

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.15031    | 0.15477    | 0.1621     |   0.6 | 81.48
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.016251   | 0.021566   | 0.02856    |   2.2 | 11.35
Output  | 0.0012321  | 0.0012329  | 0.0012347  |   0.0 |  0.65
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.01238    |            |       |  6.52

Nlocal:    50.5 ave 52 max 48 min
Histogram: 8 0 0 0 0 8 0 0 0 16
Nghost:    1360 ave 1378 max 1342 min
Histogram: 8 0 4 4 0 0 8 0 0 8
Neighs:    5629.5 ave 5800 max 5343 min
Histogram: 4 4 0 0 0 0 6 8 0 10
FullNghs:  11259 ave 11648 max 10630 min
Histogram: 8 0 0 0 0 8 0 0 0 16

Total # of neighbors = 360288
Ave neighs/atom = 222.95
Neighbor list builds = 0
Dangerous builds = 0
if "${BRVar} == 1" then "fix BoxRelax all box/relax iso 0.0 vmax 0.001" "min_style cg" "minimize 1e-08 0.0 1000 1000000" "unfix BoxRelax" "min_style cg" "minimize 1e-08 0.0 1000 1000000"
fix BoxRelax all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1e-08 0.0 1000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Per MPI rank memory allocation (min/avg/max) = 6.142 | 6.143 | 6.144 Mbytes
Step Temp PotEng c_eatoms 
     135                        0       -7146.523648687691       -7146.523648687694 
     137                        0       -7146.524967963245       -7146.524967963242 
Loop time of 0.00580474 on 32 procs for 2 steps with 1616 atoms

100.0% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -7146.52364869     -7146.52496796     -7146.52496796
  Force two-norm initial, final = 21.6883 0.143376
  Force max component initial, final = 21.6883 0.143169
  Final line search alpha, max atom move = 1 0.143169
  Iterations, force evaluations = 2 4

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.0028751  | 0.002955   | 0.0030332  |   0.1 | 50.91
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00026667 | 0.00032955 | 0.00043009 |   0.0 |  5.68
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.00252    |            |       | 43.42

Nlocal:    50.5 ave 52 max 48 min
Histogram: 8 0 0 0 0 8 0 0 0 16
Nghost:    1357.56 ave 1378 max 1342 min
Histogram: 8 0 4 4 4 0 8 0 0 4
Neighs:    5629.5 ave 5792 max 5349 min
Histogram: 4 4 0 0 0 0 6 6 2 10
FullNghs:  11259 ave 11648 max 10630 min
Histogram: 8 0 0 0 0 8 0 0 0 16

Total # of neighbors = 360288
Ave neighs/atom = 222.95
Neighbor list builds = 0
Dangerous builds = 0
unfix BoxRelax
min_style cg
minimize 1e-08 0.0 1000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Per MPI rank memory allocation (min/avg/max) = 6.142 | 6.143 | 6.144 Mbytes
Step Temp PotEng c_eatoms 
     137                        0       -7146.524967963245       -7146.524967963242 
     138                        0       -7146.524970307072       -7146.524970307071 
Loop time of 0.00316279 on 32 procs for 1 steps with 1616 atoms

98.7% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -7146.52496796     -7146.52496796     -7146.52497031
  Force two-norm initial, final = 0.0077025 0.00471959
  Force max component initial, final = 0.000611397 0.000480214
  Final line search alpha, max atom move = 1 0.000480214
  Iterations, force evaluations = 1 2

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.0016776  | 0.0017243  | 0.0017736  |   0.1 | 54.52
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00017468 | 0.00020911 | 0.00026979 |   0.0 |  6.61
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.001229   |            |       | 38.87

Nlocal:    50.5 ave 52 max 48 min
Histogram: 8 0 0 0 0 8 0 0 0 16
Nghost:    1357.56 ave 1378 max 1342 min
Histogram: 8 0 4 4 4 0 8 0 0 4
Neighs:    5629.5 ave 5792 max 5349 min
Histogram: 4 4 0 0 0 0 6 6 2 10
FullNghs:  11259 ave 11648 max 10630 min
Histogram: 8 0 0 0 0 8 0 0 0 16

Total # of neighbors = 360288
Ave neighs/atom = 222.95
Neighbor list builds = 0
Dangerous builds = 0

if "${NVTVar} == 1" then "fix 			NVTFix all nvt temp 300.0 300.0 1.0" "run 			5000" "unfix 			NVTFix" "fix BoxRelax2 all box/relax iso 0.0 vmax 0.001" "min_style cg" "minimize 1e-08 0.0 1000 1000000" "unfix BoxRelax2"
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Per MPI rank memory allocation (min/avg/max) = 5.017 | 5.018 | 5.019 Mbytes
Step Temp PotEng c_eatoms 
     138                        0       -7146.524970307071        -7146.52497030707 
Loop time of 3.46289e-06 on 32 procs for 0 steps with 1616 atoms

102.9% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 3.463e-06  |            |       |100.00

Nlocal:    50.5 ave 52 max 48 min
Histogram: 8 0 0 0 0 8 0 0 0 16
Nghost:    1357.56 ave 1378 max 1342 min
Histogram: 8 0 4 4 4 0 8 0 0 4
Neighs:    5629.5 ave 5791 max 5347 min
Histogram: 4 3 1 0 0 0 6 6 2 10
FullNghs:  11259 ave 11648 max 10630 min
Histogram: 8 0 0 0 0 8 0 0 0 16

Total # of neighbors = 360288
Ave neighs/atom = 222.95
Neighbor list builds = 0
Dangerous builds = 0
undump			GBBulkDump
shell cd ..

variable NNi equal count(nickel)
variable noNi6 equal ${NNi}
variable noNi6 equal 1608

variable NS equal count(sulfur)
variable noS6 equal ${NS}
variable noS6 equal 8

variable E equal "c_eatoms"
if "${IntPotVar} == 1" then "variable Ei equal $E" else "variable Ei equal $E/23.06035  #eV"
variable Ei equal $E
variable Ei equal -7146.52497030707

variable EIGBB equal ${Ei}
variable EIGBB equal -7146.52497030707

#---- Remove nickel to find the energy of the impurity atoms ---------
delete_atoms group nickel
Deleted 1608 atoms, new total = 8

reset_timestep	0
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
WARNING: New thermo_style command, previous thermo_modify settings will be lost (../output.cpp:708)

shell cd DumpFiles
dump MonolayerGBBulkDump all custom 5000 dump.MonolayerGBBulk id type x y z c_csym c_PE
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Per MPI rank memory allocation (min/avg/max) = 4.014 | 4.14 | 5.015 Mbytes
Step PotEng Lx Ly Lz Press Pxx Pyy Pzz c_eatoms 
       0   -1.0133422    29.008712     87.02114    7.0371609   -480.32867 -3.8095258e-24 -5.077091e-25    -1440.986   -1.0133422 
Loop time of 3.71458e-06 on 32 procs for 0 steps with 8 atoms

77.4% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 3.715e-06  |            |       |100.00

Nlocal:    0.25 ave 2 max 0 min
Histogram: 28 0 0 0 0 0 0 0 0 4
Nghost:    4.625 ave 20 max 0 min
Histogram: 24 0 0 0 0 0 0 0 2 6
Neighs:    1 ave 9 max 0 min
Histogram: 28 0 0 0 0 0 0 2 0 2
FullNghs:  2 ave 16 max 0 min
Histogram: 28 0 0 0 0 0 0 0 0 4

Total # of neighbors = 64
Ave neighs/atom = 8
Neighbor list builds = 0
Dangerous builds = 0
undump MonolayerGBBulkDump
shell cd ..

variable N equal count(all)
variable noS6Bulk equal $N
variable noS6Bulk equal 8

variable E equal "c_eatoms"
if "${IntPotVar} == 1" then "variable EiS equal $E" else "variable EiS equal $E/23.06035  #eV"
variable EiS equal $E
variable EiS equal -1.01334221788951

variable EIPGBB equal ${EiS}
variable EIPGBB equal -1.01334221788951
#---------------------------------------------


#####################################################
#################       (g)       ###################
#####################################################
#log 	log.FSBulk

# --------------- INITIALIZATION ------------------
clear
if "${IntPotVar} == 1" then "units metal" "atom_style atomic" else "units real" "atom_style charge"
units metal
atom_style atomic

dimension 3
boundary p p p

shell cd RestartFiles
read_data 	data.FSBulk
  orthogonal box = (0.004876 0.034254 0.000425) to (29.01 93.5447 7.03672)
  4 by 8 by 1 MPI processor grid
  reading atoms ...
  1616 atoms
  read_data CPU = 0.00210091 secs
shell cd ..
# ---------- Define Interatomic Potential ---------------------
if "${IntPotVar} == 1" then "pair_style      eam/alloy" else "pair_style      reax/c NULL safezone 1.6 mincap 100"
pair_style      eam/alloy

pair_coeff * * ../../${IntPot} Ni S Ni Ni
pair_coeff * * ../../NiS.eam.alloy Ni S Ni Ni

if "${IntPotVar} == 0" then "mass 1 58.71" "mass 2 32.06" "mass 1 58.71" "mass 1 58.71"

# ---------- Groups ---------------------------------
group 			nickel type 1
1608 atoms in group nickel
group 			sulfur type 2
8 atoms in group sulfur

#--------------------------- Computations ----------------------------
compute csym all centro/atom fcc
compute PE all pe/atom
compute eatoms all reduce sum c_PE

if "${IntPotVar} == 0" then "fix QEQREAX all qeq/reax 1 0.0 10.0 1.0e-6 ../../paramNiS.qeq"

# ------------ Thermo style ---------------------------------
thermo 			100
thermo_style    custom step temp pe c_eatoms
thermo_modify   format float %24.16g

# ----------- Dump Results ---------------------------------
shell cd DumpFiles
dump            FSBulkDump all custom 100 dump.FSBulk id type x y z c_PE c_csym
variable NVTVar equal 0
variable BRVar equal 1
include ../../../min.config
min_style cg
minimize 1e-25 1e-25 1000000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Neighbor list info ...
  update every 1 steps, delay 0 steps, check yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 8.5
  ghost atom cutoff = 8.5
  binsize = 4.25, bins = 7 23 2
  2 neighbor lists, perpetual/occasional/extra = 1 1 0
  (1) pair eam/alloy, perpetual
      attributes: half, newton on
      pair build: half/bin/atomonly/newton
      stencil: half/bin/3d/newton
      bin: standard
  (2) compute centro/atom, occasional
      attributes: full, newton on
      pair build: full/bin/atomonly
      stencil: full/bin/3d
      bin: standard
Per MPI rank memory allocation (min/avg/max) = 6.128 | 6.297 | 6.768 Mbytes
Step Temp PotEng c_eatoms 
       0                        0       -7102.915642150504       -7102.915642150502 
     100                        0       -7112.888897389244       -7112.888897389243 
     200                        0       -7112.888907231338       -7112.888907231335 
     215                        0        -7112.88890723135       -7112.888907231352 
Loop time of 0.3023 on 32 procs for 215 steps with 1616 atoms

100.0% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -7102.91564215     -7112.88890723     -7112.88890723
  Force two-norm initial, final = 8.93246 1.37596e-06
  Force max component initial, final = 1.9532 2.27245e-07
  Final line search alpha, max atom move = 1 2.27245e-07
  Iterations, force evaluations = 215 424

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.19576    | 0.2439     | 0.2689     |   5.2 | 80.68
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.022369   | 0.043557   | 0.089738   |  11.8 | 14.41
Output  | 0.0023673  | 0.002368   | 0.0023718  |   0.0 |  0.78
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.01247    |            |       |  4.13

Nlocal:    50.5 ave 56 max 35 min
Histogram: 4 0 4 0 0 0 0 0 0 24
Nghost:    1313.38 ave 1409 max 1069 min
Histogram: 8 0 0 0 0 0 0 0 4 20
Neighs:    5432.75 ave 6449 max 3271 min
Histogram: 8 0 0 0 0 0 0 0 16 8
FullNghs:  10865.5 ave 12544 max 6164 min
Histogram: 4 0 4 0 0 0 0 0 0 24

Total # of neighbors = 347696
Ave neighs/atom = 215.158
Neighbor list builds = 0
Dangerous builds = 0
if "${BRVar} == 1" then "fix BoxRelax all box/relax iso 0.0 vmax 0.001" "min_style cg" "minimize 1e-08 0.0 1000 1000000" "unfix BoxRelax" "min_style cg" "minimize 1e-08 0.0 1000 1000000"
fix BoxRelax all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1e-08 0.0 1000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Per MPI rank memory allocation (min/avg/max) = 6.128 | 6.297 | 6.768 Mbytes
Step Temp PotEng c_eatoms 
     215                        0       -7112.888907231349       -7112.888907231351 
     217                        0       -7112.896770646034       -7112.896770646036 
Loop time of 0.00611557 on 32 procs for 2 steps with 1616 atoms

100.0% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -7112.88890723     -7112.89677052     -7112.89677065
  Force two-norm initial, final = 52.9551 0.0731067
  Force max component initial, final = 52.9551 0.0242473
  Final line search alpha, max atom move = 0.00534169 0.000129522
  Iterations, force evaluations = 2 4

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.0024069  | 0.002999   | 0.0033088  |   0.6 | 49.04
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00026277 | 0.00053016 | 0.0010671  |   0.0 |  8.67
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.002586   |            |       | 42.29

Nlocal:    50.5 ave 56 max 35 min
Histogram: 4 0 4 0 0 0 0 0 0 24
Nghost:    1315.12 ave 1414 max 1069 min
Histogram: 8 0 0 0 0 0 0 0 4 20
Neighs:    5434.75 ave 6382 max 3283 min
Histogram: 6 2 0 0 0 0 0 0 10 14
FullNghs:  10869.5 ave 12544 max 6180 min
Histogram: 4 0 4 0 0 0 0 0 0 24

Total # of neighbors = 347824
Ave neighs/atom = 215.238
Neighbor list builds = 0
Dangerous builds = 0
unfix BoxRelax
min_style cg
minimize 1e-08 0.0 1000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Per MPI rank memory allocation (min/avg/max) = 6.128 | 6.297 | 6.768 Mbytes
Step Temp PotEng c_eatoms 
     217                        0       -7112.896770646036       -7112.896770646036 
     246                        0       -7112.902238934008       -7112.902238934009 
Loop time of 0.0420074 on 32 procs for 29 steps with 1616 atoms

100.0% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -7112.89677065     -7112.90220346     -7112.90223893
  Force two-norm initial, final = 0.0689685 0.020736
  Force max component initial, final = 0.0108528 0.00383632
  Final line search alpha, max atom move = 1 0.00383632
  Iterations, force evaluations = 29 58

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.027427   | 0.033731   | 0.03662    |   1.8 | 80.30
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0031323  | 0.0058045  | 0.011811   |   4.3 | 13.82
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.002472   |            |       |  5.88

Nlocal:    50.5 ave 56 max 35 min
Histogram: 4 0 4 0 0 0 0 0 0 24
Nghost:    1315.12 ave 1414 max 1069 min
Histogram: 8 0 0 0 0 0 0 0 4 20
Neighs:    5434.75 ave 6304 max 3336 min
Histogram: 8 0 0 0 0 0 0 0 4 20
FullNghs:  10869.5 ave 12544 max 6180 min
Histogram: 4 0 4 0 0 0 0 0 0 24

Total # of neighbors = 347824
Ave neighs/atom = 215.238
Neighbor list builds = 0
Dangerous builds = 0

if "${NVTVar} == 1" then "fix 			NVTFix all nvt temp 300.0 300.0 1.0" "run 			5000" "unfix 			NVTFix" "fix BoxRelax2 all box/relax iso 0.0 vmax 0.001" "min_style cg" "minimize 1e-08 0.0 1000 1000000" "unfix BoxRelax2"
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Per MPI rank memory allocation (min/avg/max) = 5.003 | 5.172 | 5.643 Mbytes
Step Temp PotEng c_eatoms 
     246                        0       -7112.902238934008       -7112.902238934009 
Loop time of 3.49095e-06 on 32 procs for 0 steps with 1616 atoms

123.5% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 3.491e-06  |            |       |100.00

Nlocal:    50.5 ave 56 max 37 min
Histogram: 4 0 4 0 0 0 0 0 20 4
Nghost:    1315.12 ave 1414 max 1067 min
Histogram: 8 0 0 0 0 0 0 0 4 20
Neighs:    5435.62 ave 6184 max 3341 min
Histogram: 4 4 0 0 0 0 0 0 4 20
FullNghs:  10869.5 ave 12544 max 6180 min
Histogram: 4 0 4 0 0 0 0 0 0 24

Total # of neighbors = 347824
Ave neighs/atom = 215.238
Neighbor list builds = 0
Dangerous builds = 0
undump			FSBulkDump
shell cd ..

variable NNi equal count(nickel)
variable noNi7 equal ${NNi}
variable noNi7 equal 1608

variable NS equal count(sulfur)
variable noS7 equal ${NS}
variable noS7 equal 8

variable E equal "c_eatoms"
if "${IntPotVar} == 1" then "variable Ei equal $E" else "variable Ei equal $E/23.06035  #eV"
variable Ei equal $E
variable Ei equal -7112.90223893401

variable EIFSB equal ${Ei}
variable EIFSB equal -7112.90223893401
#---- Remove nickel to find the energy of the impurity atoms ---------
delete_atoms group nickel
Deleted 1608 atoms, new total = 8

reset_timestep	0
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
WARNING: New thermo_style command, previous thermo_modify settings will be lost (../output.cpp:708)

shell cd DumpFiles
dump MonolayerFSBulkDump all custom 5000 dump.MonolayerFSBulk id type x y z c_csym c_PE
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Per MPI rank memory allocation (min/avg/max) = 4 | 4.294 | 5.018 Mbytes
Step PotEng Lx Ly Lz Press Pxx Pyy Pzz c_eatoms 
       0    -1.020065    28.996545    93.482712    7.0342094   -450.23371 -1.7129776e-10 -1.7377722e-09   -1350.7011    -1.020065 
Loop time of 3.51307e-06 on 32 procs for 0 steps with 8 atoms

129.0% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 3.513e-06  |            |       |100.00

Nlocal:    0.25 ave 2 max 0 min
Histogram: 28 0 0 0 0 0 0 0 0 4
Nghost:    4.625 ave 20 max 0 min
Histogram: 24 0 0 0 0 0 0 0 2 6
Neighs:    1 ave 9 max 0 min
Histogram: 28 0 0 0 0 0 0 2 0 2
FullNghs:  2 ave 16 max 0 min
Histogram: 28 0 0 0 0 0 0 0 0 4

Total # of neighbors = 64
Ave neighs/atom = 8
Neighbor list builds = 0
Dangerous builds = 0
undump MonolayerFSBulkDump
shell cd ..

variable N equal count(all)
variable noS7Bulk equal $N
variable noS7Bulk equal 8

variable E equal "c_eatoms"
if "${IntPotVar} == 1" then "variable EiS equal $E" else "variable EiS equal $E/23.06035  #eV"
variable EiS equal $E
variable EiS equal -1.02006504712677

variable EIPFSB equal ${EiS}
variable EIPFSB equal -1.02006504712677

#---------------------------------------------

# Write results
shell cd Results
print 			"""${noNi6} ${noS6} ${EIGBB} ${EIPGBB} ${noNi7} ${noS7} ${EIFSB} ${EIPFSB}""" append BulkEnergy screen no
shell cd ..

######################################
# SIMULATION DONE
print "All done"
All done
print "noNi6 = ${noNi6}"
noNi6 = 1608
print "noS6 = ${noS6}"
noS6 = 8
print "EIGBB = ${EIGBB} eV"
EIGBB = -7146.52497030707 eV
print "EIPGBB = ${EIPGBB} eV"
EIPGBB = -1.01334221788951 eV
print "noNi7 = ${noNi7}"
noNi7 = 1608
print "noS7 = ${noS7}"
noS7 = 8
print "EIFSB = ${EIFSB} eV"
EIFSB = -7112.90223893401 eV
print "EIPFSB = ${EIPFSB} eV"
EIPFSB = -1.02006504712677 eV
Total wall time: 0:00:00
