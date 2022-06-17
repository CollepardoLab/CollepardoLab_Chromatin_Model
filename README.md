First, compile LAMMPS with the following options:
cmake -D BUILD_SHARED_LIBS=yes -D BUILD_TOOLS=yes -D PKG_ASPHERE=yes -D PKG_RIGID=yes -D PKG_MOLECULE=yes -D PKG_PLUGIN=yes  ../cmake

At the cmake stage, must provide the path to the lammps src
cmake -DLAMMPS_HEADER_DIR=/home/rob/software/mylammps/src ../CollepardoLab_Chromatin_Model/
make install

Make sure the path to the shared library is on $LD_LIBRARY_PATH


# Source code for "Nucleosome plasticity is a critical element of chromatin liquid–liquid phase separation and multivalent nucleosome interactions"

We are delighted to share our model with the community. Please use it freely and cite our paper: https://doi.org/10.1038/s41467-021-23090-3 .
We are happy to answer any questions and comments by email (rc597@cam.ac.uk), and welcome contributions for any updates.

# System requirements

* Linux with C++ compilers (Intel and GCC tested)
* CMake 3.13 or greater
* LAMMPS (updated to Sep 29 2021) with MPI enabled
Tested on: CSD3 peta-4 cluster (https://www.hpc.cam.ac.uk/systems/peta-4) with Intel 2017 compliers and Ubuntu 21.04 with GCC 10.3.0.

# Installation guide

## To compile LAMMPS with our custom code

Download or compile a copy of LAMMPS with the following settings enabled:

* `BUILD_SHARED_LIBS=yes`
* `PKG_ASPHERE=yes`
* `PKG_RIGID=yes`
* `PKG_MOLECULE=yes`
* `PKG_PLUGIN=yes`

e.g. 

```
git clone https://github.com/lammps/lammps.git
cd <path/to/lammps>
git checkout tags/stable_29Sep2021 -b stable  
mkdir build && cd build
cmake -D BUILD_SHARED_LIBS=yes -D BUILD_TOOLS=yes -D PKG_ASPHERE=yes -D PKG_RIGID=yes -D PKG_MOLECULE=yes -D PKG_PLUGIN=yes -DCMAKE_CXX_FLAGS='-O3 -march=native -fno-math-errno -ffast-math'  ../cmake
make install
```
next, compile the LAMMPS plugin, giving it the path to the LAMMPS source folder:

```
git clone https://github.com/CollepardoLab/CollepardoLab_Chromatin_Model.git
cd <path/to/Collepardo Chromatin Model>
git checkout plugin-dev
mkdir build && cd build
cmake -DLAMMPS_HEADER_DIR=/<path/to/lammps>/src ../CollepardoLab_Chromatin_Model/
make install
```

`make install` will compile a shared object file (chromatin.so) which LAMMPS will expect to find somewhere on `$LD_LIBRARY_PATH`. Depending on where it is installed to, you may need to `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path/to/chromatin.so>`. This plugin implements the following commands in LAMMPS:

* `bond harmonic`
* `hremd_steve`
* `pair aniso`
* `pair ljlambda`

# Demo
## To run a single nucleosome system:
1. move to the "demo" directory
    
2. run with lammps
>mpirun -np 1 ./lmp -in in.run  
    
It will produce a LAMMPS trajectory file "dna.dump" this can be viewed in Ovito (https://www.ovito.org/)

Visible molecular dynamics will be observable after a few minutes runtime on a single core.
        
# Instructions to reproduce results
## To run chemically-specific 12N chromatin HREMD simulations:

1. The files are in main_simulations/input_scripts/chemically_specific_12N_165NRL_HREMD. The lammps input scripts are in.hremd_breathing and in.hremd_nonbreathing

2. run lammps using at least 16 cores
>mpirun -np 16 ./lmp -partition 16x1 -in in.hremd_breathing
    
These simulations for both breathing and non-breathing will generate the trajectories for our main results in figures 3 and 4.

## To run minimal model coexistence simulation:

1. The files are in main_simulations/input_scripts/minimal_coexistence/
    
2. to run a coexistence simulation:
>mpirun -np 16 ./lmp_DNA_mpi -in in.run

3. to reproduce the phase diagram (Figure 6) one would need to vary the parameters E1 and A in the input script. These correspond to the variables E and S respectively in table 5 of the supporting info. Additionally the breathing and non-breathing structures can be used by changing which data file is read in ("data_nonb.txt" or "data_b.txt").
