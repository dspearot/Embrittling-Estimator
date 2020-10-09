shell cd ..

include ${GBName}.config
include 100_510.config
# Interatomic Potential
if "${IntPotVar} == 1" then  "variable 		IntPot index 'NiS.eam.alloy'"  else  "variable 		IntPot index 'NiS.ff'"
variable 		IntPot index 'NiS.eam.alloy'

variable 		a0 equal 3.52
variable 		ebulkatom equal -4.44999999950988

# Lattice orientations
variable		gb1x1 equal 1
variable		gb1x2 equal -5
variable		gb1x3 equal 0
variable		gb1y1 equal 5
variable		gb1y2 equal 1
variable		gb1y3 equal 0
variable		gb1z1 equal 0
variable		gb1z2 equal 0
variable		gb1z3 equal 1
variable		gb2x1 equal -1
variable		gb2x2 equal 5
variable		gb2x3 equal 0
variable		gb2y1 equal 5
variable		gb2y2 equal 1
variable		gb2y3 equal 0
variable		gb2z1 equal 0
variable		gb2z2 equal 0
variable		gb2z3 equal -1

# Variables for S insertion
variable        deltamu equal 0.1 #eV
variable        sigma string 510
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
variable 		square1 equal 1^2+-5^2+${gb1x3}^2
variable 		square1 equal 1^2+-5^2+0^2
variable 		square2 equal ${gb1y1}^2+${gb1y2}^2+${gb1y3}^2
variable 		square2 equal 5^2+${gb1y2}^2+${gb1y3}^2
variable 		square2 equal 5^2+1^2+${gb1y3}^2
variable 		square2 equal 5^2+1^2+0^2
variable 		square3 equal ${gb1z1}^2+${gb1z2}^2+${gb1z3}^2
variable 		square3 equal 0^2+${gb1z2}^2+${gb1z3}^2
variable 		square3 equal 0^2+0^2+${gb1z3}^2
variable 		square3 equal 0^2+0^2+1^2

variable repDistx equal 2.5495097567963922
variable repDisty equal 2.5495097567963922
variable repDistz equal 1.0

# 3 Layer Structure dimensions (Step1 & Step2)
variable		xlo3Layer equal -0.001
variable		xhi3Layer equal 1*${a0}*${repDistx}
variable		xhi3Layer equal 1*3.52*${repDistx}
variable		xhi3Layer equal 1*3.52*2.54950975679639
variable		ylo3Layer equal -0.001
variable		yhi3Layer equal 8*${a0}*${repDisty}
variable		yhi3Layer equal 8*3.52*${repDisty}
variable		yhi3Layer equal 8*3.52*2.54950975679639
variable		zlo3Layer equal -0.001
variable		zhi3Layer equal 1*${a0}*${repDistz}
variable		zhi3Layer equal 1*3.52*${repDistz}
variable		zhi3Layer equal 1*3.52*1
variable        GB1 equal ${yhi3Layer}/4
variable        GB1 equal 71.7941947513863/4
variable        GB2 equal ${yhi3Layer}/4*3
variable        GB2 equal 71.7941947513863/4*3

# Iteration Numbers
variable		iter_x equal 10
variable		iter_y equal 10
variable		iter_z equal 10

variable		raster_x equal 10

# Box Dimensions for in.Structures and in.SGCMC
variable		xlo equal -0.001
variable		xhi equal 1*${a0}*${repDistx}
variable		xhi equal 1*3.52*${repDistx}
variable		xhi equal 1*3.52*2.54950975679639
variable		ylo equal -0.001
variable		yhi equal 4*${a0}*${repDisty}
variable		yhi equal 4*3.52*${repDisty}
variable		yhi equal 4*3.52*2.54950975679639
variable		zlo equal -0.001
variable		zhi equal 1*${a0}*${repDistz}
variable		zhi equal 1*3.52*${repDistz}
variable		zhi equal 1*3.52*1
variable		GB1Str equal 2*${a0}*${repDisty}
variable		GB1Str equal 2*3.52*${repDisty}
variable		GB1Str equal 2*3.52*2.54950975679639
variable		FS1Str equal 3*${a0}*${repDisty}
variable		FS1Str equal 3*3.52*${repDisty}
variable		FS1Str equal 3*3.52*2.54950975679639

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
  orthogonal box = (0.000612 0.011609 -0.000368) to (17.9449 71.7816 7.03916)
  4 by 8 by 1 MPI processor grid
  reading atoms ...
  824 atoms
  read_data CPU = 0.00331053 secs
shell cd ..
# ---------- Define Interatomic Potential ---------------------
if "${IntPotVar} == 1" then "pair_style      eam/alloy" else "pair_style      reax/c NULL safezone 1.6 mincap 100"
pair_style      eam/alloy

pair_coeff * * ../../${IntPot} Ni S Ni Ni
pair_coeff * * ../../NiS.eam.alloy Ni S Ni Ni

if "${IntPotVar} == 0" then "mass 1 58.71" "mass 2 32.06" "mass 1 58.71" "mass 1 58.71"

# ---------- Groups ---------------------------------
group 			nickel type 1
820 atoms in group nickel
group 			sulfur type 2
4 atoms in group sulfur

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
  binsize = 4.25, bins = 5 17 2
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
Per MPI rank memory allocation (min/avg/max) = 6.126 | 6.127 | 6.128 Mbytes
Step Temp PotEng c_eatoms 
       0                        0        -3640.79137385722       -3640.791373857223 
     100                        0       -3640.903104743606       -3640.903104743606 
     112                        0       -3640.903104743656       -3640.903104743656 
Loop time of 0.0978011 on 32 procs for 112 steps with 824 atoms

99.9% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -3640.79137386     -3640.90310474     -3640.90310474
  Force two-norm initial, final = 1.72282 3.18212e-06
  Force max component initial, final = 0.299716 4.38287e-07
  Final line search alpha, max atom move = 1 4.38287e-07
  Iterations, force evaluations = 112 224

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.068912   | 0.073146   | 0.077811   |   0.8 | 74.79
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.011833   | 0.016266   | 0.020577   |   1.8 | 16.63
Output  | 0.00075376 | 0.0007543  | 0.00075646 |   0.0 |  0.77
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.007634   |            |       |  7.81

Nlocal:    25.75 ave 28 max 24 min
Histogram: 8 0 0 0 0 20 0 0 0 4
Nghost:    1128.25 ave 1139 max 1114 min
Histogram: 8 0 0 0 4 0 4 4 4 8
Neighs:    2876.75 ave 3157 max 2595 min
Histogram: 4 2 2 0 2 10 8 0 2 2
FullNghs:  5753.5 ave 6190 max 5292 min
Histogram: 8 0 0 0 0 12 8 0 0 4

Total # of neighbors = 184112
Ave neighs/atom = 223.437
Neighbor list builds = 0
Dangerous builds = 0
if "${BRVar} == 1" then "fix BoxRelax all box/relax iso 0.0 vmax 0.001" "min_style cg" "minimize 1e-08 0.0 1000 1000000" "unfix BoxRelax" "min_style cg" "minimize 1e-08 0.0 1000 1000000"
fix BoxRelax all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1e-08 0.0 1000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Per MPI rank memory allocation (min/avg/max) = 6.126 | 6.127 | 6.128 Mbytes
Step Temp PotEng c_eatoms 
     112                        0       -3640.903104743656       -3640.903104743656 
     114                        0       -3640.903757639481       -3640.903757639482 
Loop time of 0.00357213 on 32 procs for 2 steps with 824 atoms

87.5% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -3640.90310474     -3640.90375758     -3640.90375764
  Force two-norm initial, final = 10.89 0.146831
  Force max component initial, final = 10.89 0.146731
  Final line search alpha, max atom move = 1 0.146731
  Iterations, force evaluations = 2 4

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.0015857  | 0.001664   | 0.0017973  |   0.1 | 46.58
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00024341 | 0.00033044 | 0.00040154 |   0.0 |  9.25
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.001578   |            |       | 44.17

Nlocal:    25.75 ave 28 max 24 min
Histogram: 8 0 0 0 0 20 0 0 0 4
Nghost:    1130.12 ave 1151 max 1114 min
Histogram: 8 0 4 0 8 2 6 0 0 4
Neighs:    2877.12 ave 3159 max 2576 min
Histogram: 2 2 2 2 2 10 8 0 2 2
FullNghs:  5752.5 ave 6190 max 5288 min
Histogram: 6 2 0 0 0 14 6 0 0 4

Total # of neighbors = 184080
Ave neighs/atom = 223.398
Neighbor list builds = 0
Dangerous builds = 0
unfix BoxRelax
min_style cg
minimize 1e-08 0.0 1000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Per MPI rank memory allocation (min/avg/max) = 6.126 | 6.127 | 6.128 Mbytes
Step Temp PotEng c_eatoms 
     114                        0       -3640.903757639481       -3640.903757639482 
     115                        0       -3640.903759140096       -3640.903759140096 
Loop time of 0.00196593 on 32 procs for 1 steps with 824 atoms

99.0% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -3640.90375764     -3640.90375764     -3640.90375914
  Force two-norm initial, final = 0.0054237 0.00417316
  Force max component initial, final = 0.000689172 0.000646462
  Final line search alpha, max atom move = 1 0.000646462
  Iterations, force evaluations = 1 2

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.00093319 | 0.00097812 | 0.0010537  |   0.0 | 49.75
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00013963 | 0.00019186 | 0.00023951 |   0.0 |  9.76
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.000796   |            |       | 40.49

Nlocal:    25.75 ave 28 max 24 min
Histogram: 8 0 0 0 0 20 0 0 0 4
Nghost:    1130.12 ave 1151 max 1114 min
Histogram: 8 0 4 0 8 2 6 0 0 4
Neighs:    2876.25 ave 3157 max 2574 min
Histogram: 2 2 2 2 0 10 10 0 2 2
FullNghs:  5752.5 ave 6190 max 5288 min
Histogram: 6 2 0 0 0 14 6 0 0 4

Total # of neighbors = 184080
Ave neighs/atom = 223.398
Neighbor list builds = 0
Dangerous builds = 0

if "${NVTVar} == 1" then "fix 			NVTFix all nvt temp 300.0 300.0 1.0" "run 			5000" "unfix 			NVTFix" "fix BoxRelax2 all box/relax iso 0.0 vmax 0.001" "min_style cg" "minimize 1e-08 0.0 1000 1000000" "unfix BoxRelax2"
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Per MPI rank memory allocation (min/avg/max) = 5.001 | 5.002 | 5.003 Mbytes
Step Temp PotEng c_eatoms 
     115                        0       -3640.903759140096       -3640.903759140096 
Loop time of 3.4718e-06 on 32 procs for 0 steps with 824 atoms

130.5% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 3.472e-06  |            |       |100.00

Nlocal:    25.75 ave 28 max 24 min
Histogram: 8 0 0 0 0 20 0 0 0 4
Nghost:    1130.12 ave 1151 max 1114 min
Histogram: 8 0 4 0 8 2 6 0 0 4
Neighs:    2876.25 ave 3157 max 2574 min
Histogram: 2 2 2 2 0 10 10 0 2 2
FullNghs:  5752.5 ave 6190 max 5288 min
Histogram: 6 2 0 0 0 14 6 0 0 4

Total # of neighbors = 184080
Ave neighs/atom = 223.398
Neighbor list builds = 0
Dangerous builds = 0
undump			GBBulkDump
shell cd ..

variable NNi equal count(nickel)
variable noNi6 equal ${NNi}
variable noNi6 equal 820

variable NS equal count(sulfur)
variable noS6 equal ${NS}
variable noS6 equal 4

variable E equal "c_eatoms"
if "${IntPotVar} == 1" then "variable Ei equal $E" else "variable Ei equal $E/23.06035  #eV"
variable Ei equal $E
variable Ei equal -3640.9037591401

variable EIGBB equal ${Ei}
variable EIGBB equal -3640.9037591401

#---- Remove nickel to find the energy of the impurity atoms ---------
delete_atoms group nickel
Deleted 820 atoms, new total = 4

reset_timestep	0
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
WARNING: New thermo_style command, previous thermo_modify settings will be lost (../output.cpp:708)

shell cd DumpFiles
dump MonolayerGBBulkDump all custom 5000 dump.MonolayerGBBulk id type x y z c_csym c_PE
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Per MPI rank memory allocation (min/avg/max) = 4 | 4.063 | 5.001 Mbytes
Step PotEng Lx Ly Lz Press Pxx Pyy Pzz c_eatoms 
       0  -0.50302666    17.946424    71.778701    7.0403827   -467.30933 -4.7717179e-25 -3.2108756e-25    -1401.928  -0.50302666 
Loop time of 3.53466e-06 on 32 procs for 0 steps with 4 atoms

37.1% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 3.535e-06  |            |       |100.00

Nlocal:    0.125 ave 2 max 0 min
Histogram: 30 0 0 0 0 0 0 0 0 2
Nghost:    4.375 ave 12 max 0 min
Histogram: 20 0 0 0 0 0 0 0 2 10
Neighs:    0.25 ave 4 max 0 min
Histogram: 30 0 0 0 0 0 0 0 0 2
FullNghs:  0.5 ave 8 max 0 min
Histogram: 30 0 0 0 0 0 0 0 0 2

Total # of neighbors = 16
Ave neighs/atom = 4
Neighbor list builds = 0
Dangerous builds = 0
undump MonolayerGBBulkDump
shell cd ..

variable N equal count(all)
variable noS6Bulk equal $N
variable noS6Bulk equal 4

variable E equal "c_eatoms"
if "${IntPotVar} == 1" then "variable EiS equal $E" else "variable EiS equal $E/23.06035  #eV"
variable EiS equal $E
variable EiS equal -0.503026663525414

variable EIPGBB equal ${EiS}
variable EIPGBB equal -0.503026663525414
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
  orthogonal box = (0.000612 0.011609 -0.000368) to (17.9449 78.2816 7.03916)
  2 by 16 by 1 MPI processor grid
  reading atoms ...
  824 atoms
  read_data CPU = 0.00151715 secs
shell cd ..
# ---------- Define Interatomic Potential ---------------------
if "${IntPotVar} == 1" then "pair_style      eam/alloy" else "pair_style      reax/c NULL safezone 1.6 mincap 100"
pair_style      eam/alloy

pair_coeff * * ../../${IntPot} Ni S Ni Ni
pair_coeff * * ../../NiS.eam.alloy Ni S Ni Ni

if "${IntPotVar} == 0" then "mass 1 58.71" "mass 2 32.06" "mass 1 58.71" "mass 1 58.71"

# ---------- Groups ---------------------------------
group 			nickel type 1
820 atoms in group nickel
group 			sulfur type 2
4 atoms in group sulfur

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
  binsize = 4.25, bins = 5 19 2
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
Per MPI rank memory allocation (min/avg/max) = 6.11 | 6.241 | 6.753 Mbytes
Step Temp PotEng c_eatoms 
       0                        0       -3614.047737249803       -3614.047737249805 
     100                        0       -3619.664879532353       -3619.664879532352 
     149                        0       -3619.664881001855       -3619.664881001856 
Loop time of 0.123194 on 32 procs for 149 steps with 824 atoms

99.9% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -3614.04773725       -3619.664881       -3619.664881
  Force two-norm initial, final = 6.63398 4.38374e-06
  Force max component initial, final = 2.00771 1.28216e-06
  Final line search alpha, max atom move = 1 1.28216e-06
  Iterations, force evaluations = 149 293

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.049549   | 0.091307   | 0.1024     |   5.0 | 74.12
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.013384   | 0.023065   | 0.061234   |   8.8 | 18.72
Output  | 0.00073898 | 0.0007394  | 0.00074146 |   0.0 |  0.60
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.008083   |            |       |  6.56

Nlocal:    25.75 ave 30 max 5 min
Histogram: 2 0 2 0 0 0 0 0 2 26
Nghost:    1051.75 ave 1182 max 778 min
Histogram: 8 0 0 0 0 0 2 2 4 16
Neighs:    2752.62 ave 3405 max 393 min
Histogram: 4 0 0 0 0 0 2 0 10 16
FullNghs:  5505.25 ave 6724 max 606 min
Histogram: 2 2 0 0 0 0 0 4 2 22

Total # of neighbors = 176168
Ave neighs/atom = 213.796
Neighbor list builds = 0
Dangerous builds = 0
if "${BRVar} == 1" then "fix BoxRelax all box/relax iso 0.0 vmax 0.001" "min_style cg" "minimize 1e-08 0.0 1000 1000000" "unfix BoxRelax" "min_style cg" "minimize 1e-08 0.0 1000 1000000"
fix BoxRelax all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1e-08 0.0 1000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Per MPI rank memory allocation (min/avg/max) = 6.11 | 6.264 | 6.753 Mbytes
Step Temp PotEng c_eatoms 
     149                        0       -3619.664881001854       -3619.664881001856 
     151                        0       -3619.672895047634       -3619.672895047633 
Loop time of 0.00370902 on 32 procs for 2 steps with 824 atoms

99.6% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
          -3619.664881      -3619.6728949     -3619.67289505
  Force two-norm initial, final = 38.1572 0.0928833
  Force max component initial, final = 38.1572 0.0521069
  Final line search alpha, max atom move = 0.00898878 0.000468377
  Iterations, force evaluations = 2 4

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.0010099  | 0.0016617  | 0.0018602  |   0.6 | 44.80
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00024254 | 0.00039738 | 0.00096908 |   0.0 | 10.71
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.00165    |            |       | 44.49

Nlocal:    25.75 ave 30 max 7 min
Histogram: 2 2 0 0 0 0 0 0 2 26
Nghost:    1052.31 ave 1163 max 780 min
Histogram: 8 0 0 0 0 0 2 2 0 20
Neighs:    2751.06 ave 3363 max 543 min
Histogram: 4 0 0 0 0 0 2 0 6 20
FullNghs:  5504.62 ave 6736 max 862 min
Histogram: 2 2 0 0 0 0 2 2 2 22

Total # of neighbors = 176148
Ave neighs/atom = 213.772
Neighbor list builds = 0
Dangerous builds = 0
unfix BoxRelax
min_style cg
minimize 1e-08 0.0 1000 1000000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:168)
Per MPI rank memory allocation (min/avg/max) = 6.11 | 6.28 | 6.753 Mbytes
Step Temp PotEng c_eatoms 
     151                        0       -3619.672895047634       -3619.672895047633 
     177                        0       -3619.678435182382       -3619.678435182381 
Loop time of 0.0227393 on 32 procs for 26 steps with 824 atoms

100.0% CPU use with 32 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -3619.67289505     -3619.67840534     -3619.67843518
  Force two-norm initial, final = 0.0768907 0.0223504
  Force max component initial, final = 0.015706 0.00517567
  Final line search alpha, max atom move = 1 0.00517567
  Iterations, force evaluations = 26 52

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.0095728  | 0.016479   | 0.018642   |   2.1 | 72.47
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0025479  | 0.0041403  | 0.010174   |   3.5 | 18.21
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.002121   |            |       |  9.33

Nlocal:    25.75 ave 30 max 7 min
Histogram: 2 2 0 0 0 0 0 0 2 26
Nghost:    1052.31 ave 1163 max 780 min
Histogram: 8 0 0 0 0 0 2 2 0 20
Neighs:    2752.31 ave 3372 max 543 min
Histogram: 4 0 0 0 0 0 2 0 6 20
FullNghs:  5504.62 ave 6736 max 862 min
Histogram: 2 2 0 0 0 0 2 2 2 22

Total # of neighbors = 176148
Ave neighs/atom = 213.772
Neighbor list builds = 0
Dangerous builds = 0

if "${NVTVar} == 1" then "fix 			NVTFix all nvt temp 300.0 300.0 1.0" "run 			5000" "unfix 			NVTFix" "fix BoxRelax2 all box/relax iso 0.0 vmax 0.001" "min_style cg" "minimize 1e-08 0.0 1000 1000000" "unfix BoxRelax2"
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Per MPI rank memory allocation (min/avg/max) = 4.985 | 5.155 | 5.628 Mbytes
Step Temp PotEng c_eatoms 
     177                        0       -3619.678435182382       -3619.678435182381 
Loop time of 3.40887e-06 on 32 procs for 0 steps with 824 atoms

124.7% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 3.409e-06  |            |       |100.00

Nlocal:    25.75 ave 30 max 7 min
Histogram: 2 2 0 0 0 0 0 0 2 26
Nghost:    1054.75 ave 1163 max 798 min
Histogram: 8 0 0 0 0 0 2 2 0 20
Neighs:    2750.94 ave 3362 max 543 min
Histogram: 4 0 0 0 0 0 2 0 6 20
FullNghs:  5504.62 ave 6736 max 862 min
Histogram: 2 2 0 0 0 0 2 2 2 22

Total # of neighbors = 176148
Ave neighs/atom = 213.772
Neighbor list builds = 0
Dangerous builds = 0
undump			FSBulkDump
shell cd ..

variable NNi equal count(nickel)
variable noNi7 equal ${NNi}
variable noNi7 equal 820

variable NS equal count(sulfur)
variable noS7 equal ${NS}
variable noS7 equal 4

variable E equal "c_eatoms"
if "${IntPotVar} == 1" then "variable Ei equal $E" else "variable Ei equal $E/23.06035  #eV"
variable Ei equal $E
variable Ei equal -3619.67843518238

variable EIFSB equal ${Ei}
variable EIFSB equal -3619.67843518238
#---- Remove nickel to find the energy of the impurity atoms ---------
delete_atoms group nickel
Deleted 820 atoms, new total = 4

reset_timestep	0
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
WARNING: New thermo_style command, previous thermo_modify settings will be lost (../output.cpp:708)

shell cd DumpFiles
dump MonolayerFSBulkDump all custom 5000 dump.MonolayerFSBulk id type x y z c_csym c_PE
run 0
WARNING: No fixes defined, atoms won't move (../verlet.cpp:52)
Per MPI rank memory allocation (min/avg/max) = 3.984 | 4.216 | 5.004 Mbytes
Step PotEng Lx Ly Lz Press Pxx Pyy Pzz c_eatoms 
       0  -0.50734449      17.9367    78.237078    7.0365681   -432.58734 -9.2483778e-21 -1.0275684e-20    -1297.762  -0.50734449 
Loop time of 3.56877e-06 on 32 procs for 0 steps with 4 atoms

113.8% CPU use with 32 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0          | 0          | 0          |   0.0 |  0.00
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0          | 0          | 0          |   0.0 |  0.00
Output  | 0          | 0          | 0          |   0.0 |  0.00
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 3.569e-06  |            |       |100.00

Nlocal:    0.125 ave 2 max 0 min
Histogram: 30 0 0 0 0 0 0 0 0 2
Nghost:    5.5 ave 18 max 0 min
Histogram: 22 0 0 0 0 0 0 0 2 8
Neighs:    0.25 ave 4 max 0 min
Histogram: 30 0 0 0 0 0 0 0 0 2
FullNghs:  0.5 ave 8 max 0 min
Histogram: 30 0 0 0 0 0 0 0 0 2

Total # of neighbors = 16
Ave neighs/atom = 4
Neighbor list builds = 0
Dangerous builds = 0
undump MonolayerFSBulkDump
shell cd ..

variable N equal count(all)
variable noS7Bulk equal $N
variable noS7Bulk equal 4

variable E equal "c_eatoms"
if "${IntPotVar} == 1" then "variable EiS equal $E" else "variable EiS equal $E/23.06035  #eV"
variable EiS equal $E
variable EiS equal -0.507344489905781

variable EIPFSB equal ${EiS}
variable EIPFSB equal -0.507344489905781

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
noNi6 = 820
print "noS6 = ${noS6}"
noS6 = 4
print "EIGBB = ${EIGBB} eV"
EIGBB = -3640.9037591401 eV
print "EIPGBB = ${EIPGBB} eV"
EIPGBB = -0.503026663525414 eV
print "noNi7 = ${noNi7}"
noNi7 = 820
print "noS7 = ${noS7}"
noS7 = 4
print "EIFSB = ${EIFSB} eV"
EIFSB = -3619.67843518238 eV
print "EIPFSB = ${EIPFSB} eV"
EIPFSB = -0.507344489905781 eV
Total wall time: 0:00:00
