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
------------------------------------------------------------------------- */

#include "temper_2_mod.h"
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
#include "RIGID/fix_rigid_small.h"
using namespace LAMMPS_NS;

// #define Temper2Mod2_DEBUG 1

/* ---------------------------------------------------------------------- */

Temper2Mod::Temper2Mod(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

Temper2Mod::~Temper2Mod()
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
   perform Temper2Moding with inter-world swaps
------------------------------------------------------------------------- */

void Temper2Mod::command(int narg, char **arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to Temper2Mod");
  if (domain->box_exist == 0)
    error->all(FLERR,"Temper2Mod command before simulation box is defined");
  if (narg != 7 && narg != 8)
    error->universe_all(FLERR,"Illegal Temper2Mod command");

  int nsteps = utils::inumeric(FLERR, arg[0], false, lmp);
  nevery = utils::inumeric(FLERR, arg[1], false, lmp);
  double temp = utils::numeric(FLERR, arg[2], false, lmp);

  // ignore Temper2Mod command, if walltime limit was already reached

  if (timer->is_timeout()) return;

  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix)
    error->universe_all(FLERR,"Temper2Moding fix ID is not defined");

 for (whichfix2 = 0; whichfix2 < modify->nfix; whichfix2++)
    if (strcmp(arg[4],modify->fix[whichfix2]->id) == 0) break;
  if (whichfix2 == modify->nfix)
    error->universe_all(FLERR,"Temper2Moding fix2 ID is not defined");
 
  for (whichfixRigid = 0; whichfixRigid < modify->nfix; whichfixRigid++)
    if (strcmp(arg[5],modify->fix[whichfixRigid]->id) == 0) break;
  if (whichfixRigid == modify->nfix)
    error->universe_all(FLERR,"Temper2Moding fix2 ID is not defined");


  seed_swap = utils::inumeric(FLERR, arg[6], false, lmp);
  seed_boltz = utils::inumeric(FLERR, arg[7], false, lmp);

  my_set_temp = universe->iworld;
  if (narg == 9) my_set_temp = utils::inumeric(FLERR, arg[8], false, lmp);
  if ((my_set_temp < 0) || (my_set_temp >= universe->nworlds))
    error->universe_one(FLERR,"Illegal Temper2Modature index");

  // swap frequency must evenly divide total # of timesteps

  if (nevery <= 0)
    error->universe_all(FLERR,"Invalid frequency in Temper2Mod command");
  nswaps = nsteps/nevery;
  if (nswaps*nevery != nsteps)
    error->universe_all(FLERR,"Non integer # of swaps in Temper2Mod command");

  // fix style must be appropriate for Temper2Modature control, i.e. it needs
  // to provide a working Fix::reset_target() and must not change the volume.

  if ((!utils::strmatch(modify->fix[whichfix]->style,"^nvt")) &&
      (!utils::strmatch(modify->fix[whichfix]->style,"^langevin")) &&
      (!utils::strmatch(modify->fix[whichfix]->style,"^gl[de]$")) &&
      (!utils::strmatch(modify->fix[whichfix]->style,"^rigid/nvt")) &&
      (!utils::strmatch(modify->fix[whichfix]->style,"^temp/")))
    error->universe_all(FLERR,"Temper2Moding Temperature fix is not supported");


 if ((!utils::strmatch(modify->fix[whichfix2]->style,"^nvt")) &&
      (!utils::strmatch(modify->fix[whichfix2]->style,"^langevin")) &&
      (!utils::strmatch(modify->fix[whichfix2]->style,"^gl[de]$")) &&
      (!utils::strmatch(modify->fix[whichfix2]->style,"^rigid/nvt")) &&
      (!utils::strmatch(modify->fix[whichfix2]->style,"^temp/")))
    error->universe_all(FLERR,"Temper2Moding Temperature fix is not supported");
 
 if ((!utils::strmatch(modify->fix[whichfixRigid]->style,"^rigid/nve/small")))
    error->universe_all(FLERR,"Temper2Moding Temperature fix is not supported");


  fix_rigid_ptr = static_cast<FixRigidSmall *> (modify->fix[whichfixRigid]);

  // setup for long Temper2Moding run

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
  if (id < 0) error->all(FLERR,"Temper2Moding could not find thermo_pe compute");
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

  // create static list of set Temper2Modatures
  // allgather Temper2Moding arg "temp" across root procs
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

  // if restarting Temper2Moding, reset temp target of Fix to current my_set_temp

  if (narg == 9) {
    double new_temp = set_temp[my_set_temp];
    modify->fix[whichfix]->reset_target(new_temp);
    modify->fix[whichfix2]->reset_target(new_temp);
  }

  // setup Temper2Moding runs

  int i,which,partner,swap,partner_set_temp,partner_world;
  double pe,pe_partner,boltz_factor,new_temp;

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up Temper2Moding ...\n");

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

  timer->init();
  timer->barrier_start();

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

    // compute PE
    // notify compute it will be called at next swap

    pe = pe_compute->compute_scalar();
    pe_compute->addstep(update->ntimestep + nevery);

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

    swap = 0;
    if (partner != -1) {
      if (me_universe > partner)
        MPI_Send(&pe,1,MPI_DOUBLE,partner,0,universe->uworld);
      else
        MPI_Recv(&pe_partner,1,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);

      if (me_universe < partner) {
        boltz_factor = (pe - pe_partner) *
          (1.0/(boltz*set_temp[my_set_temp]) -
           1.0/(boltz*set_temp[partner_set_temp]));
        if (boltz_factor >= 0.0) swap = 1;
        else if (ranboltz->uniform() < exp(boltz_factor)) swap = 1;
      }

      if (me_universe < partner)
        MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld);
      else
        MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,MPI_STATUS_IGNORE);

#ifdef Temper2Mod_DEBUG
      if (me_universe < partner)
        printf("SWAP %d & %d: yes = %d,Ts = %d %d, PEs = %g %g, Bz = %g %g\n",
               me_universe,partner,swap,my_set_temp,partner_set_temp,
               pe,pe_partner,boltz_factor,exp(boltz_factor));
#endif

    }

    // bcast swap result to other procs in my world

    MPI_Bcast(&swap,1,MPI_INT,0,world);

    // rescale kinetic energy via velocities if move is accepted

    if (swap) scale_velocities(partner_set_temp,my_set_temp);

    // if my world swapped, all procs in world reset temp target of Fix

    if (swap) {
      new_temp = set_temp[partner_set_temp];
      modify->fix[whichfix]->reset_target(new_temp);
      modify->fix[whichfix2]->reset_target(new_temp);
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

    if (me_universe == 0) print_status();
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

void Temper2Mod::scale_velocities(int t_partner, int t_me)
{
  double sfactor = sqrt(set_temp[t_partner]/set_temp[t_me]);

  double **v = atom->v;
  double **angmom = atom->angmom;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    v[i][0] = v[i][0]*sfactor;
    v[i][1] = v[i][1]*sfactor;
    v[i][2] = v[i][2]*sfactor;

    // rescale rotational
    angmom[i][0]*=sfactor;
    angmom[i][1]*=sfactor;
    angmom[i][2]*=sfactor;

  }

  // rescale (nve/small) rigid bodies
  int ibody;
  int nlocal_body = fix_rigid_ptr->nlocal_body;
  //printf("%d bodies to scale\n",nlocal_body);
  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    FixRigidSmall::Body *b = &fix_rigid_ptr->body[ibody];
    b->vcm[0]*=sfactor;
    b->vcm[1]*=sfactor;
    b->vcm[2]*=sfactor;
    
    b->conjqm[0]*=sfactor;   
    b->conjqm[1]*=sfactor;   
    b->conjqm[2]*=sfactor;   
    b->conjqm[3]*=sfactor;   
  }
}

/* ----------------------------------------------------------------------
   proc 0 prints current Temper2Moding status
------------------------------------------------------------------------- */

void Temper2Mod::print_status()
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