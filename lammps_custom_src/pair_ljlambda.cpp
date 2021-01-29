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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_ljlambda.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJLambda::PairLJLambda(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJLambda::~PairLJLambda()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    memory->destroy(lambda); //JM
  }
}

/* ---------------------------------------------------------------------- */

void PairLJLambda::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double r, rinv, screening;

  double TWO_1_3 = pow(2.0,(1.0/3.0)); //JM

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
//JM
        if (rsq < cut_coulsq[itype][jtype]) {
          r = sqrt(rsq);
          rinv = 1.0/r;
          screening = exp(-kappa*r);
          forcecoul = qqrd2e * qtmp*q[j] * screening * (kappa + rinv);
        } else forcecoul = 0.0; //JM
//        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
//      else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
//JM
    if (rsq <= TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) 
            forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          else            
            forcelj = lambda[itype][jtype] * r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]); //JM
        } else forcelj = 0.0;

        fpair = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype])
      ecoul = factor_coul * qqrd2e * qtmp*q[j] * rinv * screening;       //ecoul = factor_coul * qqrd2e * qtmp*q[j]*sqrt(r2inv);
          else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
//JM
            if (rsq <= TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype])
              evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) + (1-lambda[itype][jtype])*epsilon[itype][jtype];
            else
        evdwl = lambda[itype][jtype]*r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) - offset[itype][jtype];
//JM
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJLambda::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(cut_coul,n+1,n+1,"pair:cut_coul");
  memory->create(cut_coulsq,n+1,n+1,"pair:cut_coulsq");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(lambda,n+1,n+1,"pair:lambda"); //JM
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJLambda::settings(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR,"Illegal pair_style command");
  
  kappa = force->numeric(FLERR,arg[0]);
  cut_lj_global = force->numeric(FLERR,arg[1]);
  if (narg == 1) cut_coul_global = cut_lj_global;
  else cut_coul_global = force->numeric(FLERR,arg[2]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_lj[i][j] = cut_lj_global;
          cut_coul[i][j] = cut_coul_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJLambda::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 7) 
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);
  double lambda_one = force->numeric(FLERR,arg[4]); //JM

  double cut_lj_one = cut_lj_global;
  double cut_coul_one = cut_coul_global;
  if (narg >= 6) cut_coul_one = cut_lj_one = force->numeric(FLERR,arg[5]);
  if (narg == 7) cut_coul_one = force->numeric(FLERR,arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      lambda[i][j] = lambda_one; //JM
      cut_lj[i][j] = cut_lj_one;
      cut_coul[i][j] = cut_coul_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJLambda::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style ljlambda requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJLambda::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i],cut_lj[j][j]);
    cut_coul[i][j] = mix_distance(cut_coul[i][i],cut_coul[j][j]);
  }

  double cut = MAX(cut_lj[i][j],cut_coul[i][j]);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  lambda[j][i] = lambda[i][j]; //JM
  sigma[j][i] = sigma[i][j]; //JM
  epsilon[j][i] = epsilon[i][j]; //JM
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut_lj[i][j]*cut_lj[i][j]*cut_lj[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJLambda::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&lambda[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
        fwrite(&cut_coul[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJLambda::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&lambda[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
          fread(&cut_coul[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lambda[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_coul[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJLambda::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJLambda::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&kappa,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJLambda::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g \n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJLambda::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g \n",i,j,epsilon[i][j],sigma[i][j],cut_lj[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJLambda::single(int i, int j, int itype, int jtype,
                                double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,forcecoul,forcelj,phicoul,philj;
  double r,rinv,screening;

  double TWO_1_3 = pow(2.0,(1.0/3.0));

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq[itype][jtype]) {
    r = sqrt(rsq);
    rinv = 1.0/r;
    screening = exp(-kappa*r);
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j] * screening * (kappa + rinv);
  } else forcecoul = 0.0;
//    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j]*sqrt(r2inv);
//  else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    if (rsq <= TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype])
      forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    else 
      forcelj = lambda[itype][jtype] * r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  }
//forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
//} else forcelj = 0.0;

  fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j] * rinv * screening; //phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*sqrt(r2inv);
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    if (rsq <= TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype])
      philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) + (1-lambda[itype][jtype])*epsilon[itype][jtype];
    else
      philj = lambda[itype][jtype]*r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) - offset[itype][jtype];
//    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
//      offset[itype][jtype];
    eng += factor_lj*philj;
  }

  return eng;
}


void *PairLJLambda::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;
  return NULL;
}

/* ----------------------------------------------------------------------
   changes kappa and cutoffs for HREMD, added by Steve
------------------------------------------------------------------------- */
double PairLJLambda::change_kappa(double new_kappa)
{
  //printf("here, k = %f, setting new k to %f\n",kappa, new_kappa);
  double old_kappa = kappa;
  kappa = new_kappa;
  double new_hremd_cutoff = 3.5/new_kappa;
  
  //// call init_one() for each I,J
  //// set cutsq for each I,J, used to neighbor
  //// cutforce = max of all I,J cutoffs
  cut_coul_global = new_hremd_cutoff;
  cutforce = new_hremd_cutoff;
  double cut;

  for (int i = 1; i <= atom->ntypes; i++){
    for (int j = i; j <= atom->ntypes; j++) {
      
      if (cut_coul[i][j]>0.0000001){ // floating point zero
        cut_coul[i][j] = cut_coul[j][i] = new_hremd_cutoff;
      
        cut = init_one(i,j);
        cutsq[i][j] = cutsq[j][i] = cut*cut;
      } 
  }
  }

  return old_kappa;
}

void PairLJLambda::compute_ecoul_steve(double * old_E_tot, double * new_E_tot, double old_kappa, double new_kappa){

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double r, rinv, old_screening, new_screening;

  double TWO_1_3 = pow(2.0,(1.0/3.0)); //JM

  double old_ecoul = 0.0;
  double new_ecoul = 0.0;
  double old_hremd_cutoff = 3.5/old_kappa;
  double new_hremd_cutoff = 3.5/new_kappa;
  
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // if charged
      double qq = qtmp*q[j];
      if (fabs(qq)>0.0000001){ // floating point zero

      //if (rsq < cutsq[itype][jtype]) {
//      r2inv = 1.0/rsq;
//JM
        //if (rsq < cut_coulsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        
        if( r < old_hremd_cutoff){
           old_screening = exp(-old_kappa*r);
           old_ecoul += factor_coul * qqrd2e * qq * rinv * old_screening;       //ecoul = factor_coul * qqrd2e * qtmp*q[j]*sqrt(r2inv);
        }
        if( r < new_hremd_cutoff){
          new_screening = exp(-new_kappa*r);
          new_ecoul += factor_coul * qqrd2e * qq * rinv * new_screening;       //ecoul = factor_coul * qqrd2e * qtmp*q[j]*sqrt(r2inv);
        }
       } 
     } 
   }
  *old_E_tot = old_ecoul;
  *new_E_tot = new_ecoul;
}
