# Source code for "Nucleosome plasticity is a critical element of chromatin liquid–liquid phase separation and multivalent nucleosome interactions"

We are delighted to share our model with the community. Please use it freely and cite our paper DOI:XXX (Preprint: https://doi.org/10.1101/2020.11.23.391599 ) . 
We are happy to answer any questions and comments by email (rc597@cam.ac.uk), and welcome contributions for any updates.




# System requirements

Linux with C++ compilers with MPI. 
Tested on: CSD3 peta-4 cluster (https://www.hpc.cam.ac.uk/systems/peta-4) with Intel 2017 compliers


# Installation guide

## To compile LAMMPS with our custom code

1. clone a copy of LAMMPS
> git clone https://github.com/lammps/lammps.git

2. checkout stable version 3rd March 2020
> cd lammps  
> git checkout tags/stable_3Mar2020 -b stable  

3. copy all our code from lammps_custom_code into lammps/src

4. move Makefile_DNA_mpi from lammps/src into lammps/src/MAKE

5. Install required lammps packages
>make yes-asphere  
>make yes-rigid  
>make yes-molecule  

6. compile using our makefile, note this is for Intel compilers only
>make DNA_mpi  

7. the executable will be lmp_DNA_mpi

# Demo
## To run a single nucleosome system:
1. move to the "demo" directory
    
2. run with lammps
>mpirun -np 1 ./lmp_DNA_mpi -in in.run  
    
It will produce a LAMMPS trajectory file "dna.dump" this can be viewed in Ovito (https://www.ovito.org/)

Visible molecular dynamics will be observable after a few minutes runtime on a single core.
        
# Instructions to reproduce results
## To run chemically-specific 12N chromatin HREMD simulations:

1. The files are in main_simulations/input_scripts/chemically_specific_12N_165NRL_HREMD. The lammps input scripts are in.hremd_breathing and in.hremd_nonbreathing

2. run lammps using at least 16 cores
>mpirun -np 16 ./lmp_DNA_mpi -partition 16x1 -in in.hremd_breathing
    
These simulations for both breathing and non-breathing will generate the trajectories for our main results in figures 3 and 4.


## To run minimal model coexistence simulation:

1. The files are in main_simulations/input_scripts/minimal_coexistence/
    
2. to run a coexistence simulation:
>mpirun -np 16 ./lmp_DNA_mpi -in in.run

3. to reproduce the phase diagram (Figure 6) one would need to vary the parameters E1 and A in the input script. These correspond to the variables E and S respectively in table 5 of the supporting info. Additionally the breathing and non-breathing structures can be used by changing which data file is read in ("data_nonb.txt" or "data_b.txt").
