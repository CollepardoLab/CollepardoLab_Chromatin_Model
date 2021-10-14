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
   DNA rigid basepair interaction, modified from bond_harmonic.cpp
   S.FARR 2018
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "bond_harmonic_DNA.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "math_extra.h"
#include "math_const.h"
//#include "math_vector.h"
#include <fstream>
#include <sstream>
//#define PRINT_DEBUG
//#define NAN_DEBUG
//#define BOND_ORDER_DEBUG
//#define SLIDE_DEBUG

using namespace LAMMPS_NS;


/* ----------------------------------------------------------------------
   Functions added by S.Farr 2018
-------------------------------------------------------------------------*/

// modified guarded acos function incase it receives
// 1.0000000000000001 or similar
double BondHarmonic_DNA::arcos(double x){
  if(x>=1.0){
    x = 1.0;
  }else if(x<=-1.0){
    x=-1.0;
  }
  return acos(x);
}

double BondHarmonic_DNA::mag_vec(const double * v){
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// does v'.M.v where v' is the transpose of v
// 6d vectors (c arrays)
double BondHarmonic_DNA::vec_dot_mat_dot_vec(const double * vec, double mat[6][6]){
  double ret=0.0;

  double temp[6];
  //first do Mv:
  int i,j;
  for(i=0; i<6; i++){
    temp[i] = 0.0;
    for(j=0;j<6;j++){
      temp[i]+=mat[i][j]*vec[j];
    }
  }
  //now the dot product
  for(i=0;i<6;i++){
    ret+=vec[i]*temp[i];
  }

  return ret;
}
/// y = rotate "x" by "angle" about "axis".
/// axis needs to be unit vector
/// uses Rodrigue's Rotation Formula:
/// y = x cos(t) + (k cross x)sin(t) + k(k.x)(1-cos(t)) *
/// where t = angle, k = axis
void BondHarmonic_DNA::rotation(const double *x, double *y, double *axis, double angle){

  double kcrossx[3];
  MathExtra::cross3(axis,  x, kcrossx);
  double kdotx;
  kdotx = MathExtra::dot3(x,axis);
  double cosa = cos(angle);
  double sina = sin(angle);
  y[0] = x[0]*cosa + kcrossx[0]*sina + axis[0]*kdotx -axis[0]*kdotx*cosa;
  y[1] = x[1]*cosa + kcrossx[1]*sina + axis[1]*kdotx -axis[1]*kdotx*cosa;
  y[2] = x[2]*cosa + kcrossx[2]*sina + axis[2]*kdotx -axis[2]*kdotx*cosa;

}


/// Computes the energy of the bonds
/// @params:
///  x1 coordinates of bp1
///  q1 quaternion of bp1
///  x2 coordinates of bp2
///  q2 quaternion of bp2
///  helical_params helical_paramsstruct of the bond

double BondHarmonic_DNA::compute_bond_energy(const double * x1, double * q1, const double * x2, double * q2, struct helical_paramsstruct * helical_params){

  // turn q into direction vectors

  double ex1[3],ey1[3],ez1[3];
  double ex2[3],ey2[3],ez2[3];

  MathExtra::q_to_exyz(q1,ex1,ey1,ez1);
  MathExtra::q_to_exyz(q2,ex2,ey2,ez2);


  MathExtra::norm3(ex1);
  MathExtra::norm3(ey1);
  MathExtra::norm3(ez1);
  MathExtra::norm3(ex2);
  MathExtra::norm3(ey2);
  MathExtra::norm3(ez2);

  double mstx[3],msty[3],mstz[3];
  double Q[6]; // current helical parameters of the bond
  double dQ[6];
  
  // compute helical parameters
  compute_helical_parameters(x1,ex1,ey1,ez1,x2,ex2,ey2,ez2,Q,mstx,msty,mstz);

  // dQ
  dQ[0] = Q[0] - helical_params->means[0];
  dQ[1] = Q[1] - helical_params->means[1];
  dQ[2] = Q[2] - helical_params->means[2];
  dQ[3] = Q[3] - helical_params->means[3];
  dQ[4] = Q[4] - helical_params->means[4];
  dQ[5] = Q[5] - helical_params->means[5];

  double U = vec_dot_mat_dot_vec(dQ, helical_params->K)*0.5*E_CONV_FACTOR_STEVE;
  
  return U;

}


/// Computes DNA base pair step helical parameters
/// @params:
///  x1 coordinates of bp1
///  ex1, ey1, ez1 direction vectors of bp1
///  x2 coordinates of bp2
///  ex2, ey2, ez2 direction vectors of bp2
/// @returns: 6d vector of helical parameters
///  shift, slide, rise, tilt, roll, twist
///  the angles are in degrees
///
/// Uses SCHNAaP method
/// citation:
/// Xiang-Jun Lu, M.A. El Hassan, C.A. Hunter, Structure and conformation of helical nucleic acids: analysis program (SCHNAaP)
/// Journal of Molecular Biology, Volume 273, Issue 3, 31 October 1997, Pages 668-680, ISSN 0022-2836,
/// http://doi.org/10.1006/jmbi.1997.1346.
void BondHarmonic_DNA::compute_helical_parameters(const double * x1, const double * ex1_in, const double * ey1_in, const double * ez1_in, const double * x2, const double * ex2_in, const double * ey2_in,const  double * ez2_in, double * out,double * mstx_out, double * msty_out, double * mstz_out){



  double ex1[3] = {ex1_in[0], ex1_in[1], ex1_in[2]};
  double ex2[3] = {ex2_in[0], ex2_in[1], ex2_in[2]};
  double ey1[3] = {ey1_in[0], ey1_in[1], ey1_in[2]};
  double ey2[3] = {ey2_in[0], ey2_in[1], ey2_in[2]};
  double ez1[3] = {ez1_in[0], ez1_in[1], ez1_in[2]};
  double ez2[3] = {ez2_in[0], ez2_in[1], ez2_in[2]};





  // calculate the roll-tilt angle and axis
  // rt = ez1 x ez2
  // note if ez's are parallel this operation will not work
  double rt[3];
  double gamma;

  double mstx[3],msty[3],mstz[3];
  double omega,twist,phi,roll,tilt,Dx,Dy,Dz;
  double temp[3];

  MathExtra::cross3(ez1,ez2,rt);

  if((rt[0]*rt[0] + rt[1]*rt[1] + rt[2]*rt[2]) > 0.0){
#ifdef PRINT_DEBUG
    std::cout << "1" << std::endl;
#endif
    gamma = arcos(MathExtra::dot3(ez1,ez2)/(mag_vec(ez1)*mag_vec(ez2)));

    MathExtra::norm3(rt);

    // rotate bp1 by gamma/2 about rt, and bp2 by -gamma/2 about rt.

    rotation(ex1, ex1, rt, gamma*0.5);
    rotation(ey1, ey1, rt, gamma*0.5);
    rotation(ez1, ez1, rt, gamma*0.5);

    rotation(ex2, ex2, rt, -gamma*0.5);
    rotation(ey2, ey2, rt, -gamma*0.5);
    rotation(ez2, ez2, rt, -gamma*0.5);

    // make sure they are normalised
    MathExtra::norm3(ex1);
    MathExtra::norm3(ey1);
    MathExtra::norm3(ez1);
    MathExtra::norm3(ex2);
    MathExtra::norm3(ey2);
    MathExtra::norm3(ez2);

    //direction of the MST axis obtained by averaging and normalising the
    //rotated base pair triads
    mstx[0]=(ex1[0] + ex2[0])*0.5;
    mstx[1]=(ex1[1] + ex2[1])*0.5;
    mstx[2]=(ex1[2] + ex2[2])*0.5;


    msty[0]=(ey1[0] + ey2[0])*0.5;
    msty[1]=(ey1[1] + ey2[1])*0.5;
    msty[2]=(ey1[2] + ey2[2])*0.5;


    mstz[0]=(ez1[0] + ez2[0])*0.5;
    mstz[1]=(ez1[1] + ez2[1])*0.5;
    mstz[2]=(ez1[2] + ez2[2])*0.5;


    MathExtra::norm3(mstx);
    MathExtra::norm3(msty);
    MathExtra::norm3(mstz);

    //twist, omega, is the angle between the two transformed y axes

    omega = arcos(MathExtra::dot3(ey1, ey2)/(mag_vec(ey1)*mag_vec(ey2)));


    //sign control



    MathExtra::cross3(ey1,ey2,temp);
    if ( MathExtra::dot3(temp,mstz)<0.0) omega = -omega;

    twist = omega;

    // the angle between rt axis and msty axis is phi

    phi = arcos(MathExtra::dot3(rt,msty)/(mag_vec(rt)*mag_vec(msty)));

    //sign control
    MathExtra::cross3(rt, msty, temp);
    if( MathExtra::dot3(temp,mstz)<0.0){ phi = - phi;}

    //roll and tilt

    roll = gamma*cos(phi);
    tilt = gamma*sin(phi);

    //shift, slide, rise are the components of the relative displacement ofinline void copy3(const double *v, double *ans);
    //the two base pairs triads along the x,y,z, axes of the mst
    //Di = (r2 - r1).msti

    temp[0] = x2[0] - x1[0];
    temp[1] = x2[1] - x1[1];
    temp[2] = x2[2] - x1[2];

    Dx = MathExtra::dot3(temp, mstx);
    Dy = MathExtra::dot3(temp, msty);
    Dz = MathExtra::dot3(temp, mstz);

  }else{
#ifdef PRINT_DEBUG
    std::cout << "special" << std::endl;
#endif
    //special case
    //roll-tilt angle is zero, cannot define rolltilt axis by rt=z1xz2
    // roll and tilt are zero:
    roll = 0.0;
    tilt = 0.0;

    //twist is the angle between the two y axes and x axis
    twist = arcos(MathExtra::dot3(ex1,ex2)/(mag_vec(ex1)*mag_vec(ex2)));



    mstx[0]=(ex1[0] + ex2[0])*0.5;
    mstx[1]=(ex1[1] + ex2[1])*0.5;
    mstx[2]=(ex1[2] + ex2[2])*0.5;


    msty[0]=(ey1[0] + ey2[0])*0.5;
    msty[1]=(ey1[1] + ey2[1])*0.5;
    msty[2]=(ey1[2] + ey2[2])*0.5;


    mstz[0]=(ez1[0] + ez2[0])*0.5;
    mstz[1]=(ez1[1] + ez2[1])*0.5;
    mstz[2]=(ez1[2] + ez2[2])*0.5;


    MathExtra::norm3(mstx);
    MathExtra::norm3(msty);
    MathExtra::norm3(mstz);


    temp[0] = x2[0] - x1[0];
    temp[1] = x2[1] - x1[1];
    temp[2] = x2[2] - x1[2];

    Dx = MathExtra::dot3(temp, mstx);
    Dy = MathExtra::dot3(temp, msty);
    Dz = MathExtra::dot3(temp, mstz);
  }

  out[0] = Dx;
  out[1] = Dy;
  out[2] = Dz;
  out[3] = 180.0/MathConst::MY_PI*tilt;
  out[4] = 180.0/MathConst::MY_PI*roll;
  out[5] = 180.0/MathConst::MY_PI*twist;

#ifdef NAN_DEBUG

  for(int a=0;a<6;++a){
    if(isnan(out[a])){
      std::cout << a << std::endl;
      exit(1);
    }
  }
#endif


  mstx_out[0]=mstx[0];
  mstx_out[1]=mstx[1];
  mstx_out[2]=mstx[2];


  msty_out[0]=msty[0];
  msty_out[1]=msty[1];
  msty_out[2]=msty[2];


  mstz_out[0]=mstz[0];
  mstz_out[1]=mstz[1];
  mstz_out[2]=mstz[2];


}

/* ---------------------------------------------------------------------- */

BondHarmonic_DNA::BondHarmonic_DNA(LAMMPS *lmp) : Bond(lmp)
{
  reinitflag = 1;

  //avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
}

/* ---------------------------------------------------------------------- */

BondHarmonic_DNA::~BondHarmonic_DNA()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonic_DNA::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,rk;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  tagint *tags = atom->tag;
  //----------------------------------------------------------------------------
  // Added in by S.Farr 2018
  double **torque = atom->torque;
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double fx, fy, fz;
  double * quat1, * quat2;

  // DNA rigid base pair interaction
  // U = 1/2 dQ K dQ
  // K = 6x6 matrix
  // dQ is change of helical parameters (shift,slide,rise,tilt,roll,twist) from
  // their equilibrium values
  //
  // Method:
  // - compute dQ

  double mstx[3];
  double msty[3];
  double mstz[3];

  double Q[6]; // current helical parameters of the bond
  double dQ[6];
  //double Q_mean[6] = {0.0,0.0,3.4,0.0,0.0,30.0};
  //  double Q_mean[6] = {-0.3,       -0.3,        3.3,       -2.6,        0.3,       35.4};
  //  double K[6][6] = {{1.72017,    0.19796,    0.32533,   -0.01249,    0.00576,    0.05913},
  //                    {0.19797,    2.12618,    0.75074,   -0.00581,   -0.05309,   -0.10162},
  //                    {0.32534,    0.75074,    7.64359,   -0.18348,   -0.04547,    -0.1485},
  //                    {-0.01249,   -0.00581,   -0.18349,    0.03738,    0.00211,    0.00597},
  //                    {0.00576,   -0.05309,   -0.04547,    0.00211,    0.01961,    0.00742},
  //                    {0.05913,   -0.10162,    -0.1485,    0.00597,    0.00742,    0.02761}};
  //


  for (n = 0; n < nbondlist; n++) {


    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

#ifdef BOND_ORDER_DEBUG

    int tagi = tags[i1];
    int tagj = tags[i2];

    if(tagj < tagi and abs(tagi-tagj) < 2) {
      printf("%d %d\n", tagi, tagj);
    }
#endif
    //type = bondlist[n][2];

    //delx = x[i1][0] - x[i2][0];
    //dely = x[i1][1] - x[i2][1];
    //delz = x[i1][2] - x[i2][2];

    //rsq = delx*delx + dely*dely + delz*delz;
    //r = sqrt(rsq);
    //dr = r - r0[type];
    //rk = k[type] * dr;

    // get the quaternions
    quat1 = bonus[ellipsoid[i1]].quat;
    quat2 = bonus[ellipsoid[i2]].quat;

#ifdef NAN_DEBUG
    if(quat1 == NULL){
      std::cout << "quat1 NULL" << std::endl;
    }
    if(quat2 == NULL){
      std::cout << "quat2 NULL" << std::endl;
    }

    for(int a=0; a <4;++a){
      if(isnan(quat1[a]) || isnan(quat2[a]) ){
        std::cout << "q1 " << quat1[0] <<" "<< quat1[1] << " " << quat1[2] << " "<< quat1[3] << std::endl;
        std::cout << "q2 " << quat2[0] <<" "<< quat2[1] << " " << quat2[2] << " "<< quat2[3] << std::endl;

        exit(1);

      }
    }
#endif

    // get the basepair step name
    std::string bp1 = base_pairs[tags[i1]];
    std::string bp2 = base_pairs[tags[i2]];

    // parse the basepair step
    // bp1[0] bp2[0]
    // bp1[1] bp2[1]
    // strand 1 ..A T..
    // strand 2 ..T A..   -> AT_TA
    //
    // ..A G..
    // ..T C..  -> AG_CT
    //
    std::string name = std::string(1,bp1[0]) + bp2[0] +"_"+ bp2[1] + bp1[1];

    helical_paramsstruct bond_helical_params = helical_params_map[name];

#ifdef PRINT_DEBUG
    std::cout << "name = " << name << std::endl;
#endif
/*
#ifdef SLIDE_DEBUG
    std::cout << "name = " << name << std::endl;

    std::cout << "means:\n "  << bond_helical_params.means[0] << " " << bond_helical_params.means[1] << " " << bond_helical_params.means[2] << " " << bond_helical_params.means[3] << " " << bond_helical_params.means[4] << " " << bond_helical_params.means[5] << std::endl;
    std::cout << "K:\n";
    int i,j;
    for(i=0;i<6;i++){
        for(j=0;j<6;j++){
            std::cout << bond_helical_params.K[i][j];
        }
        std::cout<<std::endl;
     }
#endif
*/


    // turn into direction vectors

    double ex1[3],ey1[3],ez1[3];
    double ex2[3],ey2[3],ez2[3];

    MathExtra::q_to_exyz(quat1,ex1,ey1,ez1);
    MathExtra::q_to_exyz(quat2,ex2,ey2,ez2);


    MathExtra::norm3(ex1);
    MathExtra::norm3(ey1);
    MathExtra::norm3(ez1);
    MathExtra::norm3(ex2);
    MathExtra::norm3(ey2);
    MathExtra::norm3(ez2);




    // compute helical parameters
    compute_helical_parameters(x[i1],ex1,ey1,ez1,x[i2],ex2,ey2,ez2,Q,mstx,msty,mstz);


#ifdef PRINT_DEBUG
    std::cout << "q1 " << quat1[0] <<" "<< quat1[1] << " " << quat1[2] << " "<< quat1[3] << std::endl;
    std::cout << "q2 " << quat2[0] <<" "<< quat2[1] << " " << quat2[2] << " "<< quat2[3] << std::endl;
    std::cout << "Q " << Q[0] <<" "<< Q[1] << " " << Q[2] << " "<< Q[3] << " " << Q[4] << " " << Q[5] << std::endl;
#endif

    // dQ
    dQ[0] = Q[0] - bond_helical_params.means[0];
    dQ[1] = Q[1] - bond_helical_params.means[1];
    dQ[2] = Q[2] - bond_helical_params.means[2];
    dQ[3] = Q[3] - bond_helical_params.means[3];
    dQ[4] = Q[4] - bond_helical_params.means[4];
    dQ[5] = Q[5] - bond_helical_params.means[5];

    // K*dQ
    double KdQ[6];
    int ii,jj;
    for(ii=0; ii<6; ii++){
      KdQ[ii] = 0.0;
      for(jj=0;jj<6;jj++){
        KdQ[ii]+= bond_helical_params.K[ii][jj]*dQ[jj];
      }
      KdQ[ii]*=E_CONV_FACTOR_STEVE;
    }

    //f = mst*(K*dQ)

    fx = mstx[0]*KdQ[0] + msty[0]*KdQ[1] + mstz[0]*KdQ[2];
    fy = mstx[1]*KdQ[0] + msty[1]*KdQ[1] + mstz[1]*KdQ[2];
    fz = mstx[2]*KdQ[0] + msty[2]*KdQ[1] + mstz[2]*KdQ[2];



    //do it numerically by finding the energy change due to a small rotation

    // initial energy
    //U = 1/2 dQ K dQ

    //double U = vec_dot_mat_dot_vec(dQ, bond_helical_params.K)*0.5*E_CONV_FACTOR_STEVE;
    
    double U = compute_bond_energy(x[i1],quat1, x[i2], quat2, &bond_helical_params);

    //std::cout << U <<  "    " << U_new << std::endl;
    double h = 0.00001;
    double inv2h = 1.0/(2.0*h);
    double invh = 1.0/h;


    // compute the central difference 
    // T = V'(x) = [V(x+h) - V(x-h)]/2h
   


    double q_rot_x[4]  = {cos( h*0.5),sin( h*0.5),0,0};
    double q_rot_nx[4] = {cos(-h*0.5),sin(-h*0.5),0,0};
    double q_rot_y[4]  = {cos( h*0.5),0,sin( h*0.5),0};
    double q_rot_ny[4] = {cos(-h*0.5),0,sin(-h*0.5),0};
    double q_rot_z[4]  = {cos( h*0.5),0,0,sin( h*0.5)};
    double q_rot_nz[4] = {cos(-h*0.5),0,0,sin(-h*0.5)};


    MathExtra::qnormalize(q_rot_x);
    MathExtra::qnormalize(q_rot_nx);
    MathExtra::qnormalize(q_rot_y);
    MathExtra::qnormalize(q_rot_ny);
    MathExtra::qnormalize(q_rot_z);
    MathExtra::qnormalize(q_rot_nz);


    double q1_h_px[4];
    double q1_h_nx[4];
    double q1_h_py[4];
    double q1_h_ny[4];
    double q1_h_pz[4];
    double q1_h_nz[4];
    double q2_h_px[4];
    double q2_h_nx[4];
    double q2_h_py[4];
    double q2_h_ny[4];
    double q2_h_pz[4];
    double q2_h_nz[4];


    MathExtra::quatquat(q_rot_x,  quat1, q1_h_px);
    MathExtra::quatquat(q_rot_nx, quat1, q1_h_nx);
    MathExtra::quatquat(q_rot_y,  quat1, q1_h_py);
    MathExtra::quatquat(q_rot_ny, quat1, q1_h_ny);
    MathExtra::quatquat(q_rot_z,  quat1, q1_h_pz);
    MathExtra::quatquat(q_rot_nz, quat1, q1_h_nz);    
    MathExtra::quatquat(q_rot_x,  quat2, q2_h_px);
    MathExtra::quatquat(q_rot_nx, quat2, q2_h_nx);
    MathExtra::quatquat(q_rot_y,  quat2, q2_h_py);
    MathExtra::quatquat(q_rot_ny, quat2, q2_h_ny);
    MathExtra::quatquat(q_rot_z,  quat2, q2_h_pz);
    MathExtra::quatquat(q_rot_nz, quat2, q2_h_nz);



    MathExtra::qnormalize(q1_h_px);
    MathExtra::qnormalize(q1_h_nx);
    MathExtra::qnormalize(q1_h_py);
    MathExtra::qnormalize(q1_h_ny);
    MathExtra::qnormalize(q1_h_pz);
    MathExtra::qnormalize(q1_h_nz);
    MathExtra::qnormalize(q2_h_px);
    MathExtra::qnormalize(q2_h_nx);
    MathExtra::qnormalize(q2_h_py);
    MathExtra::qnormalize(q2_h_ny);
    MathExtra::qnormalize(q2_h_pz);
    MathExtra::qnormalize(q2_h_nz);



   

    double tx1 = -inv2h*(compute_bond_energy(x[i1],q1_h_px,x[i2],quat2,&bond_helical_params) - compute_bond_energy(x[i1],q1_h_nx,x[i2],quat2,&bond_helical_params));
    double ty1 = -inv2h*(compute_bond_energy(x[i1],q1_h_py,x[i2],quat2,&bond_helical_params) - compute_bond_energy(x[i1],q1_h_ny,x[i2],quat2,&bond_helical_params));
    double tz1 = -inv2h*(compute_bond_energy(x[i1],q1_h_pz,x[i2],quat2,&bond_helical_params) - compute_bond_energy(x[i1],q1_h_nz,x[i2],quat2,&bond_helical_params));

    double tx2 = -inv2h*(compute_bond_energy(x[i1],quat1,x[i2],q2_h_px,&bond_helical_params) - compute_bond_energy(x[i1],quat1,x[i2],q2_h_nx,&bond_helical_params));
    double ty2 = -inv2h*(compute_bond_energy(x[i1],quat1,x[i2],q2_h_py,&bond_helical_params) - compute_bond_energy(x[i1],quat1,x[i2],q2_h_ny,&bond_helical_params));
    double tz2 = -inv2h*(compute_bond_energy(x[i1],quat1,x[i2],q2_h_pz,&bond_helical_params) - compute_bond_energy(x[i1],quat1,x[i2],q2_h_nz,&bond_helical_params));
  


//    //printf("here\n");
//    double rot_ex1[3],rot_ey1[3],rot_ez1[3];
//    double rot_ex2[3],rot_ey2[3],rot_ez2[3];
//    double new_U;
//
//    //// compute for atom 1
//
//    // x
//    rotation(ey1,rot_ey1,ex1,h);
//    rotation(ez1,rot_ez1,ex1,h);
//    MathExtra::norm3(rot_ey1);
//    MathExtra::norm3(rot_ez1);
//
//
//    compute_helical_parameters(x[i1],ex1,rot_ey1,rot_ez1,x[i2],ex2,ey2,ez2,Q,mstx,msty,mstz);
//    dQ[0] = Q[0] - bond_helical_params.means[0];
//    dQ[1] = Q[1] - bond_helical_params.means[1];
//    dQ[2] = Q[2] - bond_helical_params.means[2];
//    dQ[3] = Q[3] - bond_helical_params.means[3];
//    dQ[4] = Q[4] - bond_helical_params.means[4];
//    dQ[5] = Q[5] - bond_helical_params.means[5];
//    new_U = vec_dot_mat_dot_vec(dQ, bond_helical_params.K)*0.5*E_CONV_FACTOR_STEVE;
//
//    double dtx1 = -(new_U - U) * invh;
//
//     // y
//    rotation(ex1,rot_ex1,ey1,h);
//    rotation(ez1,rot_ez1,ey1,h);
//    MathExtra::norm3(rot_ex1);
//    MathExtra::norm3(rot_ez1);
//
//    compute_helical_parameters(x[i1],rot_ex1,ey1,rot_ez1,x[i2],ex2,ey2,ez2,Q,mstx,msty,mstz);
//    dQ[0] = Q[0] - bond_helical_params.means[0];
//    dQ[1] = Q[1] - bond_helical_params.means[1];
//    dQ[2] = Q[2] - bond_helical_params.means[2];
//    dQ[3] = Q[3] - bond_helical_params.means[3];
//    dQ[4] = Q[4] - bond_helical_params.means[4];
//    dQ[5] = Q[5] - bond_helical_params.means[5];
//    new_U = vec_dot_mat_dot_vec(dQ, bond_helical_params.K)*0.5*E_CONV_FACTOR_STEVE;
//
//    double dty1 = -(new_U - U) * invh;
//
//    // z
//    rotation(ex1,rot_ex1,ez1,h);
//    rotation(ey1,rot_ey1,ez1,h);
//    MathExtra::norm3(rot_ex1);
//    MathExtra::norm3(rot_ey1);
//
//    compute_helical_parameters(x[i1],rot_ex1,rot_ey1,ez1,x[i2],ex2,ey2,ez2,Q,mstx,msty,mstz);
//    dQ[0] = Q[0] - bond_helical_params.means[0];
//    dQ[1] = Q[1] - bond_helical_params.means[1];
//    dQ[2] = Q[2] - bond_helical_params.means[2];
//    dQ[3] = Q[3] - bond_helical_params.means[3];
//    dQ[4] = Q[4] - bond_helical_params.means[4];
//    dQ[5] = Q[5] - bond_helical_params.means[5];
//    new_U = vec_dot_mat_dot_vec(dQ, bond_helical_params.K)*0.5*E_CONV_FACTOR_STEVE;
//
//    double dtz1 = -(new_U - U) * invh;
//
//
//    // compute for atom 2
//
//    // x
//    rotation(ey2,rot_ey2,ex2,h);
//    rotation(ez2,rot_ez2,ex2,h);
//    MathExtra::norm3(rot_ey2);
//    MathExtra::norm3(rot_ez2);
//
//    compute_helical_parameters(x[i1],ex1,ey1,ez1,x[i2],ex2,rot_ey2,rot_ez2,Q,mstx,msty,mstz);
//    dQ[0] = Q[0] - bond_helical_params.means[0];
//    dQ[1] = Q[1] - bond_helical_params.means[1];
//    dQ[2] = Q[2] - bond_helical_params.means[2];
//    dQ[3] = Q[3] - bond_helical_params.means[3];
//    dQ[4] = Q[4] - bond_helical_params.means[4];
//    dQ[5] = Q[5] - bond_helical_params.means[5];
//    new_U = vec_dot_mat_dot_vec(dQ, bond_helical_params.K)*0.5*E_CONV_FACTOR_STEVE;
//
//    double dtx2 = -(new_U - U) * invh;
//
//
//    // y
//    rotation(ex2,rot_ex2,ey2,h);
//    rotation(ez2,rot_ez2,ey2,h);
//    MathExtra::norm3(rot_ex2);
//    MathExtra::norm3(rot_ez2);
//
//    compute_helical_parameters(x[i1],ex1,ey1,ez1,x[i2],rot_ex2,ey2,rot_ez2,Q,mstx,msty,mstz);
//    dQ[0] = Q[0] - bond_helical_params.means[0];
//    dQ[1] = Q[1] - bond_helical_params.means[1];
//    dQ[2] = Q[2] - bond_helical_params.means[2];
//    dQ[3] = Q[3] - bond_helical_params.means[3];
//    dQ[4] = Q[4] - bond_helical_params.means[4];
//    dQ[5] = Q[5] - bond_helical_params.means[5];
//    new_U = vec_dot_mat_dot_vec(dQ, bond_helical_params.K)*0.5*E_CONV_FACTOR_STEVE;
//
//    double dty2 = -(new_U - U) * invh;
//
//    // z
//    rotation(ex2,rot_ex2,ez2,h);
//    rotation(ey2,rot_ey2,ez2,h);
//    MathExtra::norm3(rot_ex2);
//    MathExtra::norm3(rot_ey2);
//
//    compute_helical_parameters(x[i1],ex1,ey1,ez1,x[i2],rot_ex2,rot_ey2,ez2,Q,mstx,msty,mstz);
//    dQ[0] = Q[0] - bond_helical_params.means[0];
//    dQ[1] = Q[1] - bond_helical_params.means[1];
//    dQ[2] = Q[2] - bond_helical_params.means[2];
//    dQ[3] = Q[3] - bond_helical_params.means[3];
//    dQ[4] = Q[4] - bond_helical_params.means[4];
//    dQ[5] = Q[5] - bond_helical_params.means[5];
//    new_U = vec_dot_mat_dot_vec(dQ, bond_helical_params.K)*0.5*E_CONV_FACTOR_STEVE;
//
//    double dtz2 = -(new_U - U) * invh;
//
//
//    // We have computed DT for each dna base, this value is the dna principle frame
//    // need to convert to the space frame?
//
//    double q1[4] = {quat1[0],quat1[1],quat1[2], quat1[3]};
//    double q2[4] = {quat2[0],quat2[1],quat2[2], quat2[3]};
//
//
//    double dt1[3] = {dtx1,dty1,dtz1};
//    double dt2[3] = {dtx2,dty2,dtz2};
//    double dt1_I[3];
//    double dt2_I[3];
//
//    quat_vec_rot(dt1_I,dt1,q1);
//    quat_vec_rot(dt2_I,dt2,q2);


#ifdef PRINT_DEBUG
    std::cout << "dtx1 = "<< dtx1 << std::endl;
    std::cout << "dty1 = "<< dty1 << std::endl;
    std::cout << "dtz1 = "<< dtz1 << std::endl;

    std::cout << "dtx2 = "<< dtx2 << std::endl;
    std::cout << "dty2 = "<< dty2 << std::endl;
    std::cout << "dtz2 = "<< dtz2 << std::endl;


    std::cout << "dtx1_I = "<< dt1_I[0] << std::endl;
    std::cout << "dty1_I = "<< dt1_I[1] << std::endl;
    std::cout << "dtz1_I = "<< dt1_I[2] << std::endl;

    std::cout << "dtx2_I = "<< dt2_I[0] << std::endl;
    std::cout << "dty2_I = "<< dt2_I[1] << std::endl;
    std::cout << "dtz2_I = "<< dt2_I[2] << std::endl;

#endif
#ifdef NAN_DEBUG
    if( isnan(dt1_I[0])  || isnan(dt1_I[1])  || isnan(dt1_I[2]) ||
        isnan(dt2_I[0])  || isnan(dt2_I[1])  || isnan(dt2_I[2])){

      std::cout << "q1 " << quat1[0] <<" "<< quat1[1] << " " << quat1[2] << " "<< quat1[3] << std::endl;
      std::cout << "q2 " << quat2[0] <<" "<< quat2[1] << " " << quat2[2] << " "<< quat2[3] << std::endl;


      std::cout << "dtx1 = "<< dt1_I[0] << std::endl;
      std::cout << "dty1 = "<< dt1_I[1] << std::endl;
      std::cout << "dtz1 = "<< dt1_I[2] << std::endl;

      std::cout << "dtx2 = "<< dt2_I[0] << std::endl;
      std::cout << "dty2 = "<< dt2_I[1] << std::endl;
      std::cout << "dtz2 = "<< dt2_I[2] << std::endl;
    }
#endif




    // force & energy

    fbond = sqrt(fx*fx+fy*fy+fz*fz);

    if (eflag) ebond = U;

    // apply force to each of 2 atoms


    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fx;
      f[i1][1] += fy;
      f[i1][2] += fz;
      torque[i1][0] += tx1;
      torque[i1][1] += ty1;
      torque[i1][2] += tz1;


    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= fx;
      f[i2][1] -= fy;
      f[i2][2] -= fz;
      torque[i2][0] += tx2;
      torque[i2][1] += ty2;
      torque[i2][2] += tz2;
    }

//    if (newton_bond || i1 < nlocal) {
//      f[i1][0] += delx*fbond;
//      f[i1][1] += dely*fbond;
//      f[i1][2] += delz*fbond;
//    }
//
//    if (newton_bond || i2 < nlocal) {
//      f[i2][0] -= delx*fbond;
//      f[i2][1] -= dely*fbond;
//      f[i2][2] -= delz*fbond;
//    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,fx,fy,fz);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonic_DNA::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondHarmonic_DNA::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  //force->bounds(FLERR,arg[0],atom->nbondtypes,ilo,ihi); // deprecated
  utils::bounds(FLERR, arg[0], 0, atom->nbondtypes, ilo, ihi, error);

  //double k_one = force->numeric(FLERR,arg[1]); // deprecated
  //double r0_one = force->numeric(FLERR,arg[2]); // deprecated

  double k_one = utils::numeric(FLERR, 0, arg[1], lmp);
  double r0_one = utils::numeric(FLERR, 0, arg[2], lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}


void BondHarmonic_DNA::init_style() {
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR, "Atom style ellipsoid required");

//  // read in the data we need
  get_helical_param_map("NAFlex_params.txt");
  get_basepairs("DNA_sequence.txt");
}

void BondHarmonic_DNA::get_helical_param_map(std::string fname){
  std::ifstream infile;
  infile.open(fname);
  if (!infile) {
    std::cerr << "Unable to open file " << fname << std::endl;
    exit(1);   // call system to stop
  }


  std::string line;

  int l = 1; //keep track of line number

  std::cout << "Reading: " << fname << std::endl;


  //title line
  if (getline(infile, line)) {
    l++;
  } else {
    std::cerr << "Error getting line " << l << " of " << fname << std::endl;
    exit(1);
  }

  // loop over the parameters
  int k = 0;
  while (getline(infile, line)) {
    l++; //empty line
    helical_paramsstruct read_params;
    std::string type;

    //get the base-pair step type
    if (getline(infile, line)) {
      l++;
      //std::cout << line << std::endl;
    } else {
      std::cerr << "Error getting line " << l << " of " << fname << std::endl;
      exit(1);
    }


    type = line;

    // read in the means
    if (getline(infile, line)) {
      l++;
      //std::cout << line << std::endl;
    } else {
      std::cerr << "Error getting line " << l << " of " << fname << std::endl;
      exit(1);
    }
    std::istringstream iss(line);
    iss >> read_params.means[0] >> read_params.means[1] >> read_params.means[2] >> read_params.means[3]
        >> read_params.means[4] >> read_params.means[5];

    //std::cout << read_params.means[0] << read_params.means[1] << read_params.means[2] << read_params.means[3] << read_params.means[4] << read_params.means[5] <<std::endl;


    // read empty line
    if (getline(infile, line)) {
      l++;
    } else {
      std::cerr << "Error getting line " << l << " of " << fname << std::endl;
      exit(1);
    }

    // read in the stiffness matrix
    for (int i = 0; i < 6; ++i) {
      if (getline(infile, line)) {
        l++;
        //std::cout << line << std::endl;
      } else {
        std::cerr << "Error getting line " << l << " of " << fname << std::endl;
        exit(1);
      }
      std::istringstream iss1(line);
      iss1 >> read_params.K[i][0] >> read_params.K[i][1] >> read_params
              .K[i][2] >> read_params.K[i][3] >> read_params.K[i][4] >> read_params.K[i][5];
      //std::cout << read_params.K[i][0] << read_params.K[i][1] << read_params.K[i][2] << read_params.K[i][3] << read_params.K[i][4] << read_params.K[i][5] << std::endl;

    }
    ++k;

    //put into map
    helical_params_map.insert({type, read_params});

  }
  if (k < 16) {
    std::cerr << "Read in less than the expected 16 helical parameter types\nProgram may not run correctly"
              << std::endl;
  }

  infile.close();

}
void BondHarmonic_DNA::get_basepairs(std::string fname ){

  std::ifstream infile;
  infile.open(fname);

  if (!infile) {
    std::cerr << "Unable to open file " << fname << std::endl;
    exit(1);   // call system to stop
  }

  std::cout << "Reading: " << fname << std::endl;

  std::string line;
  std::getline(infile, line);

  std::istringstream line1(line);

  char temp;
  int num_bps;
  line1 >> temp >> num_bps;
  for (int i = 0; i < num_bps; ++i) {
    std::string new_line;
    std::getline(infile, new_line);
    std::istringstream newline1(new_line);
    tagint id;
    std::string bp;
    newline1 >> id >> bp;
    base_pairs.insert({id, bp});
  }

  infile.close();

}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondHarmonic_DNA::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondHarmonic_DNA::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondHarmonic_DNA::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondHarmonic_DNA::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],r0[i]);
}

/* ---------------------------------------------------------------------- */

double BondHarmonic_DNA::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k[type] * dr;
  fforce = 0;
  if (r > 0.0) fforce = -2.0*rk/r;
  std::cout << "here" << std::endl;
  exit(1);
  return rk*dr;
}

/* ----------------------------------------------------------------------
    Return ptr to internal members upon request.
------------------------------------------------------------------------ */
void *BondHarmonic_DNA::extract( char *str, int &dim )
{
  dim = 1;
  if( strcmp(str,"kappa")==0) return (void*) k;
  if( strcmp(str,"r0")==0) return (void*) r0;
  return NULL;
}


