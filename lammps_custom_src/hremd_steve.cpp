/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mark Sears (SNL)

   Modifed from temper.cpp by S FARR (CAMBRIDGE) in 2019
------------------------------------------------------------------------- */

#include "hremd_steve.h"
#include <cmath>
#include <cstring>
#include "universe.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "fix.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "error.h"
#include "utils.h"
#include "pair.h"
#include "pair_ljlambda.h"
#include "neighbor.h"
using namespace LAMMPS_NS;

//#define TEMPER_DEBUG

/* ---------------------------------------------------------------------- */

HremdSteve::HremdSteve(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

HremdSteve::~HremdSteve()
{
  MPI_Comm_free(&roots);
  if (ranswap) delete ranswap;
  delete ranboltz;
  delete [] set_temp;
  delete [] temp2world;
  delete [] world2temp;
  delete [] world2root;
}

/* ----------------------------------------------------------------------
   perform tempering with inter-world swaps
------------------------------------------------------------------------- */

void HremdSteve::command(int narg, char **arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to temper");
  if (domain->box_exist == 0)
    error->all(FLERR,"HremdSteve command before simulation box is defined");
  if (narg != 7 && narg != 8)
    error->universe_all(FLERR,"Illegal temper command");

  int nsteps = force->inumeric(FLERR,arg[0]);
  nevery = force->inumeric(FLERR,arg[1]);
  double temp = force->numeric(FLERR,arg[2]);

  // ignore temper command, if walltime limit was already reached

  if (timer->is_timeout()) return;

  /// remove these bits, dont need any more
  // for (whichfix = 0; whichfix < modify->nfix; whichfix++)
  //   if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
  // if (whichfix == modify->nfix)
  //   error->universe_all(FLERR,"HremdSteveing fix ID is not defined");
  
  // for (whichfix2 = 0; whichfix2 < modify->nfix; whichfix2++)
  //   if (strcmp(arg[4],modify->fix[whichfix2]->id) == 0) break;
  // if (whichfix2 == modify->nfix)
  //   error->universe_all(FLERR,"HremdSteveing fix ID is not defined");


  seed_swap = force->inumeric(FLERR,arg[5]);
  seed_boltz = force->inumeric(FLERR,arg[6]);

  my_set_temp = universe->iworld;
  if (narg == 8) my_set_temp = force->inumeric(FLERR,arg[7]);
  if ((my_set_temp < 0) || (my_set_temp >= universe->nworlds))
    error->universe_one(FLERR,"Illegal temperature index");

  // swap frequency must evenly divide total # of timesteps

  if (nevery <= 0)
    error->universe_all(FLERR,"Invalid frequency in temper command");
  nswaps = nsteps/nevery;
  if (nswaps*nevery != nsteps)
    error->universe_all(FLERR,"Non integer # of swaps in temper command");

  // fix style must be appropriate for temperature control, i.e. it needs
  // to provide a working Fix::reset_target() and must not change the volume.

  // if ((!utils::strmatch(modify->fix[whichfix]->style,"^nvt")) &&
  //     (!utils::strmatch(modify->fix[whichfix]->style,"^langevin")) &&
  //     (!utils::strmatch(modify->fix[whichfix]->style,"^gl[de]$")) &&
  //     (!utils::strmatch(modify->fix[whichfix]->style,"^rigid/nvt")) &&
  //     (!utils::strmatch(modify->fix[whichfix]->style,"^temp/")))
  //   error->universe_all(FLERR,"HremdSteveing temperature fix is not supported");

  // setup for long tempering run

  update->whichflag = 1;
  timer->init_timeout();

  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // local storage

  me_universe = universe->me;
  MPI_Comm_rank(world,&me);
  nworlds = universe->nworlds;
  iworld = universe->iworld;
  boltz = force->boltz;

  // pe_compute = ptr to thermo_pe compute
  // notify compute it will be called at first swap

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"HremdSteveing could not find thermo_pe compute");
  Compute *pe_compute = modify->compute[id];
  pe_compute->addstep(update->ntimestep + nevery);

  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(universe->uworld,color,0,&roots);

  // RNGs for swaps and Boltzmann test
  // warm up Boltzmann RNG

  if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
  else ranswap = NULL;
  ranboltz = new RanPark(lmp,seed_boltz + me_universe);
  for (int i = 0; i < 100; i++) ranboltz->uniform();

  // world2root[i] = global proc that is root proc of world i

  world2root = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world);

  // create static list of set temperatures
  // allgather tempering arg "temp" across root procs
  // bcast from each root to other procs in world

  set_temp = new double[nworlds];
  if (me == 0) MPI_Allgather(&temp,1,MPI_DOUBLE,set_temp,1,MPI_DOUBLE,roots);
  MPI_Bcast(set_temp,nworlds,MPI_DOUBLE,0,world);

  // create world2temp only on root procs from my_set_temp
  // create temp2world on root procs from world2temp,
  //   then bcast to all procs within world

  world2temp = new int[nworlds];
  temp2world = new int[nworlds];
  if (me == 0) {
    MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
    for (int i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
  }
  MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);

  

  // setup tempering runs

  
  int i,which,partner,swap,partner_set_temp,partner_world;
  double pe,pe_partner,boltz_factor,new_temp;

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up tempering ...\n");

  printf("kappa on %d is %f\n",me_universe,set_temp[my_set_temp]);

  update->integrate->setup(1);

  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->uscreen," T%d",i);
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->ulogfile," T%d",i);
      fprintf(universe->ulogfile,"\n");
    }
    print_status();
  }
  
  PairLJLambda * mypair = dynamic_cast<PairLJLambda*>(force->pair);
  // if restarting tempering, reset temp target of Fix to current my_set_temp
  //if (narg == 7) {
  //  double new_temp = set_temp[my_set_temp];
  //  modify->fix[whichfix]->reset_target(new_temp);
  //}
  //
  // setup if we are restarting
  if (narg == 8){
    double new_temp = set_temp[my_set_temp];
    mypair->change_kappa(new_temp);
    lmp->neighbor->hremd_modify_cutoff(3.5/new_temp);
    lmp->neighbor->force_rebuild_hremd=true;
  }
 
  timer->init();
  timer->barrier_start();

  //bool did_i_swap = false; 
  //mypair->change_kappa(set_temp[my_set_temp]);
  for (int iswap = 0; iswap < nswaps; iswap++) {

    // run for nevery timesteps

    timer->init_timeout();



    update->integrate->run(nevery);
    

    // check for timeout across all procs

    int my_timeout=0;
    int any_timeout=0;
    if (timer->is_timeout()) my_timeout=1;
    MPI_Allreduce(&my_timeout, &any_timeout, 1, MPI_INT, MPI_SUM, universe->uworld);
    if (any_timeout) {
      timer->force_timeout();
      break;
    }


    // which = which of 2 kinds of swaps to do (0,1)

    if (!ranswap) which = iswap % 2;
    else if (ranswap->uniform() < 0.5) which = 0;
    else which = 1;

    // partner_set_temp = which set temp I am partnering with for this swap

    if (which == 0) {
      if (my_set_temp % 2 == 0) partner_set_temp = my_set_temp + 1;
      else partner_set_temp = my_set_temp - 1;
    } else {
      if (my_set_temp % 2 == 1) partner_set_temp = my_set_temp + 1;
      else partner_set_temp = my_set_temp - 1;
    }

    // partner = proc ID to swap with
    // if partner = -1, then I am not a proc that swaps

    partner = -1;
    if (me == 0 && partner_set_temp >= 0 && partner_set_temp < nworlds) {
      partner_world = temp2world[partner_set_temp];
      partner = world2root[partner_world];
    }

    // swap with a partner, only root procs in each world participate
    // hi proc sends PE to low proc
    // lo proc make Boltzmann decision on whether to swap
    // lo proc communicates decision back to hi proc
    /*
    double old_pe = pe_compute->compute_scalar();
    
    // change H
    PairLJLambda * mypair = dynamic_cast<PairLJLambda*>(force->pair);
    double og_kappa;
    update->integrate->cleanup();
    if (partner != -1){
      og_kappa = mypair->change_kappa(set_temp[partner_set_temp]);
    }

    // run for 0
    update->integrate->setup(1);
    update->integrate->run(0);
    
    // compute PE
    // notify compute it will be called at next swap
    
    pe = pe_compute->compute_scalar();
    pe_compute->addstep(update->ntimestep + nevery);

    double deltaPe = pe - old_pe;
    if (partner!=-1){
      printf("delta pe = %f\n",deltaPe);
    }
    pe=deltaPe;

    */
    // compute Es
    double my_old_e_coul=0.0;
    double my_new_e_coul=0.0;

    mypair->compute_ecoul_steve(&my_old_e_coul,&my_new_e_coul,set_temp[my_set_temp],set_temp[partner_set_temp]);

    double old_e_coul=0.0;
    double new_e_coul=0.0;

    MPI_Allreduce(&my_old_e_coul,&old_e_coul,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&my_new_e_coul,&new_e_coul,1,MPI_DOUBLE,MPI_SUM,world);

   
    #ifdef TEMPER_DEBUG
      printf("my E = %f\n",old_e_coul);
    #endif
    pe=new_e_coul-old_e_coul;


    // TODO: try and swap
    swap = 0;
    if (partner != -1) {
      if (me_universe > partner)
        MPI_Send(&pe,1,MPI_DOUBLE,partner,0,universe->uworld);
      else
        MPI_Recv(&pe_partner,1,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);

      if (me_universe < partner) {
          // boltz_factor = (pe - pe_partner) *
          //(1.0/(boltz*set_temp[my_set_temp]) -
          // 1.0/(boltz*set_temp[partner_set_temp]));
          boltz_factor = -(pe+pe_partner)/(boltz*300.0);
        if (boltz_factor >= 0.0) swap = 1;
        else if (ranboltz->uniform() < exp(boltz_factor)) swap = 1;
        //printf("prob = %f\n", exp(boltz_factor));
      }

      if (me_universe < partner)
        MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld);
      else
        MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,MPI_STATUS_IGNORE);

#ifdef TEMPER_DEBUG
      if (me_universe < partner)
        printf("SWAP %d & %d: yes = %d,Ts = %d %d, PEs = %g %g, Bz = %g %g\n",
               me_universe,partner,swap,my_set_temp,partner_set_temp,
               pe,pe_partner,boltz_factor,exp(boltz_factor));
#endif

    }

    // bcast swap result to other procs in my world

    MPI_Bcast(&swap,1,MPI_INT,0,world);


    //update->integrate->cleanup();
    //if (!swap && partner != -1){
    //  mypair->change_kappa(og_kappa);
    //}
    //update->integrate->setup(1);
    
    // rescale kinetic energy via velocities if move is accepted

    //if (swap) scale_velocities(partner_set_temp,my_set_temp);

    // if my world swapped, all procs in world reset temp target of Fix

    if (swap) {
      new_temp = set_temp[partner_set_temp];
      //modify->fix[whichfix]->reset_target(new_temp);
      mypair->change_kappa(new_temp);
      lmp->neighbor->hremd_modify_cutoff(3.5/new_temp);
      lmp->neighbor->force_rebuild_hremd=true;
      //reset_velocities(300.0);
      //did_i_swap = true;
    }

    // update my_set_temp and temp2world on every proc
    // root procs update their value if swap took place
    // allgather across root procs
    // bcast within my world

    if (swap) my_set_temp = partner_set_temp;
    if (me == 0) {
      MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
      for (i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
    }
    MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);

    // print out current swap status

    //printf("kappa on %d is %f\n",me_universe,set_temp[my_set_temp]);
    if (me_universe == 0) print_status();
    
    #ifdef TEMPER_DEBUG
      if(me==0) printf("world %d, cutnmax = %f\n",me_universe, lmp->neighbor->cutneighmax);
    #endif
  }

  timer->barrier_stop();

  update->integrate->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   scale kinetic energy via velocities a la Sugita
------------------------------------------------------------------------- */

void HremdSteve::scale_velocities(int t_partner, int t_me)
{
  double sfactor = sqrt(set_temp[t_partner]/set_temp[t_me]);

  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    v[i][0] = v[i][0]*sfactor;
    v[i][1] = v[i][1]*sfactor;
    v[i][2] = v[i][2]*sfactor;
  }
}



//void HremdSteve::reset_velocities(double T)
//{
//
//  double **v = atom->v;
//  double *m = atom->mass;
//  int nlocal = atom->nlocal;
//
//  for (int i = 0; i < nlocal; i++) {
//    v[i][0] = sqrt(boltz*T/m[i]) * ranboltz->gaussian();
//    v[i][1] = sqrt(boltz*T/m[i]) * ranboltz->gaussian();
//    v[i][2] = sqrt(boltz*T/m[i]) * ranboltz->gaussian();
//  }
//}
//




/* ----------------------------------------------------------------------
   proc 0 prints current tempering status
------------------------------------------------------------------------- */

void HremdSteve::print_status()
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->uscreen," %d",world2temp[i]);
    fprintf(universe->uscreen,"\n");
  }
  if (universe->ulogfile) {
    fprintf(universe->ulogfile,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->ulogfile," %d",world2temp[i]);
    fprintf(universe->ulogfile,"\n");
    fflush(universe->ulogfile);
  }
}
